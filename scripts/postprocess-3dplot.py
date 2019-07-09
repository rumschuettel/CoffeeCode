#!/usr/bin/env python

import numpy as np
import numexpr as ne
import tarfile, gzip, json
import re, math
import os, sys, multiprocessing


LOG2INV = 1./np.log(2.0)
FLOAT_TYPE = np.float64

# calculate shannon entropy numexpr
def shannon_entropy_numexpr_from_lambda(lambda_list, kSys, fac, qs=["(1.-q1-q2-q3)","q1","q2","q3"]):
    entropy = ""
    global_mult = 0

    for mult, poly in lambda_list:
        if len(poly) == 0:
            continue

        term = []
        for [coeff, [e1, e2, e3]] in poly:
            term.append(f"{fac} * {coeff} * ({qs[0]}**{(kSys - e1 - e2 - e3)}) * ({qs[1]}**{e1}) * ({qs[2]}**{e2}) * ({qs[3]}**{e3})")

        term = "(" + "+".join(term) + ")"

        # numexpr optimizes this call to remove double-evaluation
        entropy += f"-{mult} * where({term}==0,0,{term}*log({term})*{LOG2INV})"

        # accumulate multiplicity
        global_mult += mult

    return entropy, global_mult



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
    parser.add_argument("--min", metavar="MIN", type=float, default=0.0)
    parser.add_argument("--max", metavar="MAX", type=float, default=0.25)
    parser.add_argument("--step", metavar="STEP", type=float, default=2**-5)
    parser.add_argument(
        "infile",
        metavar="INFILE",
        type=str,
        help="tar archive with CoffeCode multinomial outputs in .json.gz format",
    )

    args = parser.parse_args()

    THREADS = args.threads
    assert THREADS > 0, "invalid thread count"

    PATH_THISFILE = os.path.dirname(os.path.realpath(__file__))
    INFILE = os.path.join(PATH_THISFILE, args.infile)
    assert os.path.exists(INFILE), "INFILE does not exist"

    MIN = args.min
    MAX = args.max
    STEP = args.step
    assert MIN < MAX and (MAX - MIN) / STEP > 0.0, "invalid MIN, MAX or STEP argument"

    # extract kSys from filename
    kSys_matches = re.findall(r"-kSys(\d+)-", args.infile)
    assert (
        len(kSys_matches) == 1
    ), "invalid filename, must encode kSys in the form -kSys3-"
    KSYS = int(kSys_matches[0])
    print("KSYS", KSYS)

    # build q lookup table to clone
    # this is not meant to be memory-efficient, otherwise we'd create them on the fly
    def qs_ref():
        for q1 in np.arange(MIN, MAX+STEP, STEP):
            for q2 in np.arange(MIN, MAX+STEP, STEP):
                for q3 in np.arange(MIN, MAX+STEP, STEP):
                    if q1+q2+q3 > 1:
                        continue

                    # antidegradability condition
                    if 2*((1-q1-q2)**2 + (q1-q2)**2 + (q1+q2)**2 + (1-q1-q2-2*q3)**2) - 16*math.sqrt(q1*q2*q3*(1-q1-q2-q3)) < 2:
                        continue

                    yield [q1, q2, q3]

    q1,q2,q3 = np.transpose(np.array(list(qs_ref())))

    # best CI
    best_ci = np.zeros(q1.shape, dtype=FLOAT_TYPE)
    best_ci.fill(-10**9)
    # best graph, we store the adjacency matrix
    best_graph = np.zeros(q1.shape, dtype=np.uint16)


    # iterate over files in INFILE
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
            ADJM = int(adjm_matches[0], 2)

            # only .gz.json files left
            print("processing", file_meta.name, f"with kSys={KSYS}, kEnv={KENV}")

            # read content as json
            with archive.extractfile(file_meta) as file:
                content = json.load(gzip.GzipFile(fileobj=file))

            S_expr, mult = shannon_entropy_numexpr_from_lambda(
                content["lambda"], KSYS, 1.0
            )
            S_expr_a, mult_a = shannon_entropy_numexpr_from_lambda(
                content["lambda_a"], KSYS, 2.0 ** -KENV
            )

            S = ne.evaluate(S_expr)
            S_a = ne.evaluate(S_expr_a)
            ci = (S_a - S) * (1./np.log2(mult_a))

            # update best ci and best graph
            best_graph[best_ci < ci] = ADJM
            best_ci = np.maximum(best_ci, ci)

    # output best graph and best ci
    OUTFILE = (
        os.path.join(PATH_THISFILE, args.outfile) if args.outfile else f"{INFILE}-best.npz"
    )
    print(f"saving best ci and graphs in {OUTFILE}")
    np.savez_compressed(OUTFILE, params=[MIN, MAX, STEP], ci=best_ci, graph=best_graph)

    print("done")
