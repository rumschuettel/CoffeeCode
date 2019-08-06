#!/usr/bin/env python

import numpy as np
import tarfile, gzip, json
import re, math
import os

from scipy.special import xlogy

LOG2E = math.log2(math.e)

def xlog2x(term):
    return xlogy(term, term) * LOG2E

import mpmath as mp
tompf = np.frompyfunc(mp.mpf, 1, 1)

def fsum(iterable):
    acc = None
    for a in iterable:
        if acc is None:
            acc = tompf(a)
            continue
        acc += tompf(a)
    return acc.astype(np.float64)


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
def shannon_entropy_from_lambda(lambda_list, kSys, fac, qs):
    entropy = np.zeros(qs.shape[:-1])

    def lambda_iterator():
        for mult, poly in lambda_list:
            if len(poly) == 0:
                continue

            def poly_iterator():
                for [coeff, [e1, e2, e3]] in poly:
                    yield (
                        fac
                        * coeff
                        * np.prod(qs ** (e1, e2, e3), -1)
                        * (1 - np.sum(qs, -1)) ** (kSys - e1 - e2 - e3)
                    )

            term = fsum(poly_iterator())

            # filter out zeros for x log x
            yield -mult * xlog2x(term)

    entropy = fsum(lambda_iterator())

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
        ret_a_gt0_w = ret_a_gt0.repeat(3).reshape(-1, 3)
        ret_b_lt0_w = ret_b_lt0.repeat(3).reshape(-1, 3)
        return {
            "qs": np.where(
                ret_a_gt0_w, np.where(ret_b_lt0_w, (qs_a + qs_b) / 2, qs_b), qs_a
            ),
            "val": np.where(
                ret_a_gt0, np.where(ret_b_lt0, (val_a + val_b) / 2, val_b), val_a
            ),
        }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="CoffeeCode Postprocess: MapReduce CI over grid with max reductor"
    )
    parser.add_argument(
        "--outfile",
        metavar="OUTFILE",
        type=str,
        default="",
        help="outfile prefix (default: INFILE-best.npz)",
    )
    parser.add_argument("--precision", metavar="PRECISION", type=int, default=mp.mp.dps)
    parser.add_argument("--radius", metavar="RADIUS", type=float, default=.50)
    parser.add_argument("--resolution", metavar="RESOLUTION", type=int, default=100)
    parser.add_argument("--bisections", metavar="BISECTIONS", type=int, default=20)
    parser.add_argument(
        "infile",
        metavar="INFILE",
        type=str,
        help="tar archive with CoffeCode multinomial outputs in .json.gz format",
    )
    parser.add_argument("--kTotMax", metavar="KTOTMAX", type=int, default=10000)

    args = parser.parse_args()

    PATH_THISFILE = os.path.dirname(os.path.realpath(__file__))
    INFILE = os.path.join(PATH_THISFILE, args.infile)
    assert os.path.exists(INFILE), "INFILE does not exist"

    PRECISION = args.precision
    assert PRECISION > 1, "wrong precision given"
    mp.mp.dps = PRECISION

    RADIUS = args.radius
    RESOLUTION = args.resolution
    assert RADIUS > 0 and RESOLUTION > 0, "radius and resolution need to be positive"

    BISECTIONS = args.bisections
    assert (
        BISECTIONS >= 0 and BISECTIONS < 64
    ), "number of bisections has to be between 0 and 63"

    KTOTMAX = args.kTotMax
    assert KTOTMAX > 0, "kTotMax needs to be a postive number"

    # extract kSys from filename
    kSys_matches = re.findall(r"-kSys(\d+)-", args.infile)
    if len(kSys_matches) == 1:
        KSYS = int(kSys_matches[0])
        print("archive with adjacency matrices detected. KSYS", KSYS)
    else:
        a_in_b_matches = re.findall(r"\.(\d+)-in-(\d+)\.", args.infile)
        assert (
            len(a_in_b_matches) == 1
        ), "invalid filename, must encode kSys in the form -kSys3- or .a-in-b."
        KSYS = None
        print("archive with catcodes detected. KSYS unset")

    # build q table on a spherical surface
    # todo: find a better spaced method
    def qs_radial():
        φ = np.arange(0, math.pi / 2, math.pi / 2 / RESOLUTION)
        θ = np.arange(0, math.pi / 2, math.pi / 2 / RESOLUTION)

        cos_φ = np.cos(φ)
        sin_φ = np.sin(φ)
        cos_θ = np.cos(θ)
        sin_θ = np.sin(θ)

        q1 = RADIUS * np.outer(sin_θ, cos_φ).flatten()
        q2 = RADIUS * np.outer(sin_θ, sin_φ).flatten()
        q3 = RADIUS * np.outer(cos_θ, np.ones_like(φ)).flatten()

        return np.stack((q1, q2, q3)).transpose()

    qs = qs_radial()

    # best CI
    best_ci = np.zeros(qs.shape[:-1])
    best_ci -= 10 ** 2
    # best qs
    best_qs = np.zeros_like(qs)
    # best graph, we store the adjacency matrix
    best_graph = np.zeros(best_ci.shape, dtype=np.int64)

    # iterate over files in INFILE
    adjm_ctr = 0
    with tarfile.open(INFILE, mode="r") as archive:
        for file_meta in archive:
            if file_meta.type != tarfile.REGTYPE:
                continue

            # extract kEnv
            adjm_matches = re.findall(r"([0-1]+)\.json\.gz", file_meta.name)
            if len(adjm_matches) == 1:
                KTOT = math.sqrt(len(adjm_matches[0]))
                assert KTOT.is_integer(), "adjacency matrix not square"
                KTOT = int(KTOT)

            else:
                a_in_b_matches = re.findall(r"\.(\d+)-in-(\d+)\.gz", file_meta.name)
                assert (
                    len(a_in_b_matches) == 1
                ), "invalid archived name, must encode adjm in 01 format as in 0110.json.gz or a catcode with filename ending in .a-in-b"
                # for catcodes we generally assume that KENV=1
                KSYS = int(a_in_b_matches[0][0]) * int(a_in_b_matches[0][1])
                KTOT = KSYS + 1

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

            def ci_fun(qqs):
                S = shannon_entropy_from_lambda(
                    content["lambda"], KSYS, 1.0, qqs
                )
                S_a = shannon_entropy_from_lambda(
                    content["lambda_a"], KSYS, 2.0 ** -KENV, qqs
                )
                mult_a = multiplicity_from_lambda(content["lambda_a"])

                return (S_a - S) * (1. / math.log2(mult_a))

            res = find_zero_passing_bisect(
                f=ci_fun, depth=BISECTIONS, qs_a=np.zeros_like(qs), qs_b=qs.copy()
            )

            # update best ci and best graph
            adjm_ctr += 1
            mask = np.sum(res["qs"], -1) > np.sum(best_qs, -1)
            best_graph[mask] = adjm_ctr
            best_ci[mask] = res["val"][mask]
            best_qs[mask] = res["qs"][mask]

            print("updated", mask.sum(dtype=np.int64))

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
        ci=best_ci.astype(np.float64),
        qs=best_qs.astype(np.float64),
        graph=best_graph.astype(np.float64),
    )

    print("done")

