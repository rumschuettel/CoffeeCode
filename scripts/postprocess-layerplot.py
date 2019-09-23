#!/usr/bin/env python

import numpy as np
import tarfile, gzip, json
import re, math
import sys, os, multiprocessing
import time

import torch

if torch.cuda.is_available():
    DEVICE = torch.device("cuda")
else:
    DEVICE = torch.device("cpu")

print("running on device", DEVICE)

FLOAT_TYPE = torch.float64
FLOAT_TYPE_NP = np.float64

torch.set_default_tensor_type(torch.DoubleTensor)

# get total multiplicity from lambda vectors
def multiplicity_from_lambda(lambda_list):
    global_mult = 0

    for mult, poly in lambda_list:
        if len(poly) == 0:
            continue

        # accumulate multiplicity
        global_mult += mult

    return global_mult


# calculate shannon entropy from a lambda vector
# fac can be used to rescale the terms for lambda_a
# qs is a np array with last dimension 4
def shannon_entropy_from_lambda(lambda_list, kSys, fac, qs, stable_summer=False):
    def lambda_iterator():
        for mult, poly in lambda_list:
            if len(poly) == 0:
                continue

            def poly_iterator():
                for [coeff, [e1, e2, e3]] in poly:
                    yield (
                        fac
                        * coeff
                        * torch.prod(
                            qs
                            ** torch.tensor(
                                [e1, e2, e3], dtype=FLOAT_TYPE, device=DEVICE
                            ),
                            -1,
                        )
                        * (1 - torch.sum(qs, -1)) ** (kSys - e1 - e2 - e3)
                    )

            term = sum(poly_iterator())
            # filter out zeros for x log x
            yield -mult * torch.where(term != 0, term * torch.log2(term), term)

    entropy = sum(lambda_iterator())
    return entropy


if __name__ == "__main__":
    import argparse
    import subprocess

    parser = argparse.ArgumentParser(
        description="CoffeeCode Postprocess: Calculate CI rate layers"
    )
    parser.add_argument(
        "--threads", metavar="THREADS", type=int, default=multiprocessing.cpu_count()
    )
    parser.add_argument(
        "--outfile",
        metavar="OUTFILE",
        type=str,
        default="",
        help="outfile prefix (default: INFILE-rate-layers.npz)",
    )
    parser.add_argument("--radius", metavar="RADIUS", type=float, default=0.50)
    parser.add_argument("--resolution", metavar="RESOLUTION", type=int, default=512)
    parser.add_argument("--layers", metavar="LAYERS", type=int, default=16)
    parser.add_argument("--external", metavar="EXTERNAL", type=bool, default=False)
    parser.add_argument(
        "infile",
        metavar="INFILE",
        type=str,
        help="tar archive with CoffeCode multinomial outputs in .json.gz format",
    )
    parser.add_argument("--kTotMax", metavar="KTOTMAX", type=int, default=10000)

    args = parser.parse_args()

    THREADS = args.threads
    assert THREADS > 0, "invalid thread count"
    torch.set_num_threads(THREADS)

    PATH_THISFILE = os.path.dirname(os.path.realpath(__file__))
    INFILE = os.path.join(PATH_THISFILE, args.infile)
    assert os.path.exists(INFILE), "INFILE does not exist"

    RADIUS = args.radius
    RESOLUTION = args.resolution
    LAYERS = args.layers
    assert RADIUS > 0 and RESOLUTION > 0 and LAYERS > 0, "radius, resolution and layers need to be positive"

    KTOTMAX = args.kTotMax
    assert KTOTMAX > 0, "kTotMax needs to be a postive number"

    EXTERNAL = args.external

    # extract kSys from filename
    ARCHIVE_FORMAT = None
    kTot_matches = re.findall(r"-kTot(\d+)", args.infile)
    kSys_matches = re.findall(r"-kSys(\d+)-", args.infile)
    if len(kSys_matches) == 1 and len(kTot_matches) == 0:
        KSYS = int(kSys_matches[0])
        print("archive with adjacency matrices detected. KSYS", KSYS)
        ARCHIVE_FORMAT = "adjm"

    elif len(kSys_matches) == 1 and len(kTot_matches) == 1:
        KSYS = int(kSys_matches[0])
        KTOT = int(kTot_matches[0])
        print("archive with KSYS", KSYS, "KTOT", KTOT)
        ARCHIVE_FORMAT = "full"

    else:
        a_in_b_matches = re.findall(r"\.(\d+)-in-(\d+)\.", args.infile)
        assert (
            len(a_in_b_matches) == 1
        ), "invalid filename, must encode kSys in the form -kSys3- or .a-in-b."
        print("archive with catcodes detected. KSYS unset")
        ARCHIVE_FORMAT = "cat"

    # print some info
    if EXTERNAL:
        print("calling external program pp3d or pp3dw")

    # build q table for layers
    # we build 
    def qs_layers():
        nx = np.linspace(0, RADIUS, RESOLUTION)
        ny = np.linspace(0, RADIUS, RESOLUTION)
        nz = np.linspace(0, RADIUS, LAYERS)
        
        q1, q2, q3 = np.meshgrid(nx, ny, nz)       

        # from coord x nx x ny x LAYER
        # to   LAYER x nx x ny x coord
        qs = np.stack((q1, q2, q3)).transpose((3, 1, 2, 0)).reshape(-1, 3)

        # of the 3d layers, we mask out everything that is antidegrabable
        # using the precise formula makes almost no difference in the number of points
        # so we use a simple heuristic and mask out all but the .5 unit simplex
        qs_mask = np.array([ True if (x+y+z < .5) else False for (x, y, z) in qs ])

        return qs, qs_mask

    # from the origin
    qs_full, qs_mask = qs_layers()
    qs = torch.from_numpy(qs_full[qs_mask]).contiguous().to(device=DEVICE)

    print("full qs grid:", len(qs_full), "of which not masked:", len(qs))

    # best CI
    best_ci = torch.zeros(qs.shape[:-1], dtype=FLOAT_TYPE, device=DEVICE)
    best_ci -= 10 ** 2
    # best graph, we store the adjacency matrix
    best_graph = torch.zeros(best_ci.shape, dtype=torch.int32, device=DEVICE)

    # iterate over files in INFILE
    adjm_ctr = 0
    with tarfile.open(INFILE, mode="r") as archive:
        for file_meta in archive:
            if file_meta.type != tarfile.REGTYPE:
                continue

            # extract kEnv
            if ARCHIVE_FORMAT == "adjm":
                adjm_matches = re.findall(r"([0-1]+)\.json\.gz", file_meta.name)
                assert len(adjm_matches) == 1, "cannot find adjacency matrix in filename as in 0110.json.gz"
                KTOT = math.sqrt(len(adjm_matches[0]))
                assert KTOT.is_integer(), "adjacency matrix not square"
                KTOT = int(KTOT)

            elif ARCHIVE_FORMAT == "cat":
                a_in_b_matches = re.findall(r"\.(\d+)-in-(\d+)\.json\.gz", file_meta.name)
                assert (
                    len(a_in_b_matches) == 1
                ), "cannot find catcode signature with filename ending in .a-in-b.json"
                # for catcodes we generally assume that KENV=1
                KSYS = int(a_in_b_matches[0][0]) * int(a_in_b_matches[0][1])
                KTOT = KSYS + 1
            
            else: # ARCHIVE_FORMAT == "full"
                pass

            KENV = KTOT - KSYS
            assert KENV > 0 and KENV <= KSYS, "0 < kEnv <= kSys not satisfied"

            if KTOT > KTOTMAX:
                # print("skipping", file_meta.name, f"with kSys={KSYS}, kEnv={KENV}")
                continue

            # only .gz.json files left
            print("processing", file_meta.name, f"with kSys={KSYS}, kEnv={KENV}")

            # read content as json
            with archive.extractfile(file_meta) as file:
                content = json.load(gzip.GzipFile(fileobj=file))

            multiplicity = multiplicity_from_lambda(content["lambda_a"])

            # internal ci_fun
            def ci_fun(qqs):
                print("ci_fun running... ", end="")
                sys.stdout.flush()
                start = time.time()
                S = shannon_entropy_from_lambda(
                    content["lambda"], KSYS, 1.0, qqs
                )
                S_a = shannon_entropy_from_lambda(
                    content["lambda_a"], KSYS, 2.0 ** -KENV, qqs
                )
                print(f"done; T={round(time.time() - start)} sec")

                return (S_a - S) / np.log2(multiplicity)

            # external ci_fun
            def ci_fun_ext(qqs):
                json_in = {
                    "lambda": content["lambda"],
                    "lambda_a": content["lambda_a"],
                    "qs": qqs.cpu().tolist()
                }

                process = subprocess.Popen(
                    ["pp3d/pp3d" if KTOT <= 32 else "pp3d/pp3dw", str(KSYS), str(KENV)],
                    stdout=subprocess.PIPE,
                    stdin=subprocess.PIPE,
                    encoding="ascii"
                )
                (output, err) = process.communicate(input=json.dumps(json_in))
                assert err == None, "external command failure"

                try:
                    json_out = json.loads(output)["S_a-S"]
                except:
                    print(output)
                    return

                return torch.tensor(json_out, dtype=FLOAT_TYPE) / np.log2(multiplicity)

            ci = (ci_fun if not EXTERNAL else ci_fun_ext)(qs)

            # update best ci and best graph
            adjm_ctr += 1
            mask = ci > best_ci
            best_graph[mask] = adjm_ctr
            best_ci[mask] = ci

            print("updated", mask.sum(dtype=torch.int64))

    # output best graph and best ci
    OUTFILE = (
        os.path.join(PATH_THISFILE, args.outfile)
        if args.outfile
        else f"{INFILE}-rate-layers.npz"
    )
    print(f"saving best ci, qs and graphs in {OUTFILE}")

    # first we reshape back to having a regular grid for best_ci and best_graph
    best_ci_full = torch.zeros(qs_full.shape[:-1], dtype=best_ci.dtype, device=DEVICE)
    best_ci_full[qs_mask] = best_ci
    best_ci_full[~qs_mask] = -100.

    best_graph_full = torch.zeros(qs_full.shape[:-1], dtype=best_graph.dtype, device=DEVICE)
    best_graph_full[qs_mask] = best_graph
    best_graph_full[~qs_mask] = -100

    np.savez_compressed(
        OUTFILE,
        params=[RADIUS, RESOLUTION, LAYERS],
        qs=qs_full.reshape(LAYERS, RESOLUTION, RESOLUTION, -1),
        ci=best_ci_full.reshape(LAYERS, RESOLUTION, RESOLUTION).cpu().numpy(),
        graph=best_graph_full.reshape(LAYERS, RESOLUTION, RESOLUTION).cpu().numpy(),
    )

    print("done")
