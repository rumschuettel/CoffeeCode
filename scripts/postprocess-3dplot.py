#!/usr/bin/env python

import numpy as np
import tarfile, gzip, json
import re, math
import os, multiprocessing

import torch

FLOAT_TYPE = torch.float64
FLOAT_TYPE_NP = np.float64


# get total multiplicity from lambda vectors
def multiplicity_from_lambda(lambda_list):
    global_mult = 0.

    for mult, poly in lambda_list:
        if len(poly) == 0:
            continue

        # accumulate multiplicity
        global_mult += mult

    return global_mult


# calculate shannon entropy from a lambda vector
# fac can be used to rescale the terms for lambda_a
# qs is a np array with last dimension 4
def shannon_entropy_from_lambda(lambda_list, kSys, fac, qs):
    entropy = torch.zeros(qs.shape[:-1], dtype=FLOAT_TYPE)

    for mult, poly in lambda_list:
        if len(poly) == 0:
            continue

        term = torch.zeros_like(entropy)
        for [coeff, [e1, e2, e3]] in poly:
            term += (
                fac
                * coeff
                * torch.prod(qs ** torch.tensor([e1, e2, e3], dtype=FLOAT_TYPE), -1)
                * (1 - torch.sum(qs, -1)) ** (kSys - e1 - e2 - e3)
            )

        # filter out zeros for x log x
        entropy -= mult * torch.where(term != 0, term * torch.log2(term), term)

    return entropy


# bisect q vectors to find zero passing in interval
def find_zero_passing_bisect(f, depth, qs_a, qs_b, val_a=None, val_b=None):
    if val_a is None or val_b is None:
        val_a = f(qs_a)
        val_b = f(qs_b)

    qs_midpoint = (qs_a + qs_b) / 2
    val_midpoint = f(qs_midpoint)

    mask = val_midpoint < 0
    qs_a[1 - mask] = qs_midpoint[1 - mask]
    qs_b[mask] = qs_midpoint[mask]
    val_a[1 - mask] = val_midpoint[1 - mask]
    val_b[mask] = val_midpoint[mask]

    if depth > 0:
        return find_zero_passing_bisect(f, depth - 1, qs_a, qs_b, val_a, val_b)
    else:
        # either return upper, lower, or midpoint qs and val
        ret_a_gt0 = val_a > 0
        ret_b_lt0 = val_b < 0
        ret_a_gt0_w = ret_a_gt0.repeat(3).reshape(3, -1).transpose(0, 1)
        ret_b_lt0_w = ret_b_lt0.repeat(3).reshape(3, -1).transpose(0, 1)
        return {
            "qs": torch.where(
                ret_a_gt0_w, torch.where(ret_b_lt0_w, (qs_a + qs_b) / 2, qs_b), qs_a
            ),
            "val": torch.where(
                ret_a_gt0, torch.where(ret_b_lt0, (val_a + val_b) / 2, val_b), val_a
            ),
        }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="CoffeeCode Postprocess: MapReduce CI over grid with max reductor"
    )
    parser.add_argument(
        "--threads", metavar="THREADS", type=int, default=multiprocessing.cpu_count()
    )
    parser.add_argument(
        "--outfile",
        metavar="OUTFILE",
        type=str,
        default="",
        help="outfile prefix (default: INFILE-best.npz)",
    )
    parser.add_argument("--radius", metavar="RADIUS", type=float, default=.50)
    parser.add_argument("--resolution", metavar="RESOLUTION", type=int, default=100)
    parser.add_argument("--bisections", metavar="BISECTIONS", type=int, default=20)
    parser.add_argument(
        "infile",
        metavar="INFILE",
        type=str,
        help="tar archive with CoffeCode multinomial outputs in .json.gz format",
    )

    args = parser.parse_args()

    THREADS = args.threads
    assert THREADS > 0, "invalid thread count"
    torch.set_num_threads(THREADS)

    PATH_THISFILE = os.path.dirname(os.path.realpath(__file__))
    INFILE = os.path.join(PATH_THISFILE, args.infile)
    assert os.path.exists(INFILE), "INFILE does not exist"

    RADIUS = args.radius
    RESOLUTION = args.resolution
    assert RADIUS > 0 and RESOLUTION > 0, "radius and resolution need to be positive"

    BISECTIONS = args.bisections
    assert (
        BISECTIONS > 0 and BISECTIONS < 64
    ), "number of bisections has to be between 1 and 63"

    # extract kSys from filename
    kSys_matches = re.findall(r"-kSys(\d+)-", args.infile)
    assert (
        len(kSys_matches) == 1
    ), "invalid filename, must encode kSys in the form -kSys3-"
    KSYS = int(kSys_matches[0])
    print("KSYS", KSYS)

    # build q table on a spherical surface
    # todo: find a better spaced method
    def qs_radial():
        φ = np.arange(0, math.pi / 2, math.pi / 2 / RESOLUTION, dtype=FLOAT_TYPE_NP)
        θ = np.arange(0, math.pi / 2, math.pi / 2 / RESOLUTION, dtype=FLOAT_TYPE_NP)

        cos_φ = np.cos(φ)
        sin_φ = np.sin(φ)
        cos_θ = np.cos(θ)
        sin_θ = np.sin(θ)

        q1 = RADIUS * np.outer(sin_θ, cos_φ).flatten()
        q2 = RADIUS * np.outer(sin_θ, sin_φ).flatten()
        q3 = RADIUS * np.outer(cos_θ, np.ones_like(φ)).flatten()

        return np.stack((q1, q2, q3)).transpose()

    qs = torch.from_numpy(qs_radial()).contiguous()

    # best CI
    best_ci = torch.zeros(qs.shape[:-1], dtype=FLOAT_TYPE)
    best_ci -= 10 ** 9
    # best qs
    best_qs = torch.zeros_like(qs)
    # best graph, we store the adjacency matrix
    best_graph = torch.zeros(best_ci.shape, dtype=torch.int32)

    # iterate over files in INFILE
    adjm_ctr = 0
    with tarfile.open(INFILE, mode="r") as archive:
        for file_meta in archive:
            if file_meta.type != tarfile.REGTYPE:
                continue

            # extract kEnv
            adjm_matches = re.findall(r"([0-1]+)\.json\.gz", file_meta.name)
            assert (
                len(adjm_matches) == 1
            ), "invalid archived name, must encode adjm in 01 format as in 0110.json.gz"
            KTOT = math.sqrt(len(adjm_matches[0]))
            assert KTOT.is_integer(), "adjacency matrix not square"
            KTOT = int(KTOT)
            KENV = KTOT - KSYS
            assert KENV > 0 and KENV <= KSYS, "0 < kEnv <= kSys not satisfied"

            # only .gz.json files left
            print("processing", file_meta.name, f"with kSys={KSYS}, kEnv={KENV}")

            # read content as json
            with archive.extractfile(file_meta) as file:
                content = json.load(gzip.GzipFile(fileobj=file))

            def ci_fun(qqs):
                S = shannon_entropy_from_lambda(content["lambda"], KSYS, 1.0, qqs)
                S_a = shannon_entropy_from_lambda(
                    content["lambda_a"], KSYS, 2.0 ** -KENV, qqs
                )
                mult_a = multiplicity_from_lambda(content["lambda_a"])

                return (S_a - S) * (1. / np.log2(mult_a))

            res = find_zero_passing_bisect(
                f=ci_fun, depth=BISECTIONS, qs_a=torch.zeros_like(qs), qs_b=qs
            )

            # update best ci and best graph
            adjm_ctr += 1
            mask = torch.sum(res["qs"], -1) > torch.sum(best_qs, -1)
            best_graph[mask] = adjm_ctr
            best_ci[mask] = res["val"][mask]
            best_qs[mask] = res["qs"][mask]

    # output best graph and best ci
    OUTFILE = (
        os.path.join(PATH_THISFILE, args.outfile)
        if args.outfile
        else f"{INFILE}-best.npz"
    )
    print(f"saving best ci, qs and graphs in {OUTFILE}")
    np.savez_compressed(
        OUTFILE,
        params=[RADIUS, RESOLUTION],
        ci=best_ci.numpy(),
        qs=best_qs.numpy(),
        graph=best_graph.numpy(),
    )

    print("done")
