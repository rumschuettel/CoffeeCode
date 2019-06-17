#!/usr/bin/env python3
# needs networkx, igraph
# as well as geng and vcolg from nauty27

import numpy as np 
from scipy.stats import rankdata

def compress_vertices(vertices, offset):
    """
    compresses a list of vertices
    0 1 2 3 --> 0 1 2 3
    1 1 1 1 --> 0 0 0 0
    5 7 6 6 --> 0 2 1 1
    luckily scipy.stats.rankdata with method "dense" does just that
    """
    return (rankdata(vertices, "dense").astype("int") - 1 + offset).tolist()

def compress_2nd_vertices(edges, offset):
    "applies compress_vertices to the second element of the tuples in the list"
    v1s, v2s = zip(*edges)
    v2s_compressed = compress_vertices(list(v2s), offset)
    return list(map(tuple, zip(v1s, v2s_compressed)))


## GRAPH
import networkx
import igraph
import base64

def to_graph6(graph, with_header=False):
    "igraph to graph6 ascii-encoded bytes"
    edges = graph.get_edgelist()
    nx_graph = networkx.Graph(edges)
    return networkx.readwrite.to_graph6_bytes(nx_graph, header=with_header).strip()

def from_graph6(graph6_str):
    "graph6 ascii-encoded bytes to igraph"
    nx_graph = networkx.readwrite.from_graph6_bytes(b">>graph6<<" + graph6_str.strip().encode("ascii"))
    edges = list(nx_graph.edges)
    return igraph.Graph(edges)

def canonicalize(graph):
    "canonical graph image under BLISS"
    return graph.permute_vertices(graph.canonical_permutation())

def canonicalize_short(graph):
    "short byte representation of edge list of canonical graph"
    return bytes(sum(set(canonicalize(graph).get_edgelist()), ()))

def graph_to_adj(graph):
    "return adjacency matrix as a flattened string of 0s and 1s"
    return "".join([ str(i) for r in graph.get_adjacency() for i in r ])

def filename_for_graph6(graph6_str):
    return base64.b32encode(graph6_str.encode("ascii")).decode("ascii")

import itertools
def powerset(iterable):
    "powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1, len(s)+1))


def all_possible_edges(kSys, kEnv):
    "iterator over all possible edge combinations between kSys system and kEnv environment vertices"
    vsys = range(0, kSys)
    venv = range(kSys, kSys+kEnv)
    return itertools.product(vsys, venv)

def all_possible_environment_links(kSys, kEnv):
    "iterator over all possible system-environment edge lists"
    for edges in powerset(all_possible_edges(kSys, kEnv)):
        yield list(edges)


## IO
import subprocess
import os

PATH_THISFILE = os.path.dirname(os.path.realpath(__file__))
PATH_GENG = os.path.join(PATH_THISFILE, "bin/geng")
PATH_VCOLG = os.path.join(PATH_THISFILE, "bin/vcolg")

def pipe(cmd, *args, encoding="ascii", inputline=None):    
    proc = subprocess.Popen(
        [cmd, *[str(a) for a in args]],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE if inputline != None else None
    )

    if inputline != None:
        proc.stdin.write((inputline + "\n").encode(encoding))
        proc.stdin.close()

    while True:
        line = proc.stdout.readline()
        if not line:
            break
        
        yield line.decode(encoding).strip()

    proc.wait()


## POOR MAN'S FILE SYSTEM LOCKS
PATH_FILESYSTEM_LOCK = os.path.join(PATH_THISFILE, "locks/")

def filesystem_indicator_fn(graph6_str, extension):
    return os.path.join(PATH_FILESYSTEM_LOCK, filename_for_graph6(graph6_str) + "." + extension)

def filesystem_lock(graph6_str):
    "check whether we already processed graph by touching a file indicator"
    lockfile = filesystem_indicator_fn(graph6_str, "lock")
    if os.path.isfile(lockfile):
        return False

    open(lockfile, "w").close()
    return True

def filesystem_unlock(graph6_str):
    "removes lock file indicator"
    lockfile = filesystem_indicator_fn(graph6_str, "lock")
    try:
        os.remove(lockfile)
        return True
    except:
        return False

def filesystem_is_done(graph6_str):
    "checks if that graph was processed"
    donefile = filesystem_indicator_fn(graph6_str, "done")
    return os.path.isfile(donefile)

def filesystem_mark_done(graph6_str):
    "marks a graph as done"
    donefile = filesystem_indicator_fn(graph6_str, "done")
    open(donefile, "w").close()



## NAUTY LINKS
def generate_all_graphs(vertex_count):
    for line in pipe(PATH_GENG, "-c", "-q", vertex_count):
        yield line, from_graph6(line)




## PROGRAM
import multiprocessing.pool
from functools import partial


if __name__ == "__main__":
    kSys = 3
    max_kEnv = kSys

    # preprocess environment links
    filtered_links = [
        (
            # list of edges like [(0, 3), (1, 3), (1, 4), (1, 5)]
            list(l),

            # maximum environment vertex
            max(l, key=lambda x: x[1])[1]

        ) for l in set([
            tuple(compress_2nd_vertices(edges, kSys)) for edges in all_possible_environment_links(kSys, max_kEnv)
        ])
    ]


    def get_all_nonisomorphic_environments(graph):
        "collected in a list since we will call this in a subthread"
        seen = set()
        candidates = []

        for links, max_env_vertex in filtered_links:
            kEnv = max_env_vertex+1 - kSys
            # create copy of graph with environment vertices added
            copy = graph.copy()
            copy.add_vertices(kEnv)
            copy.add_edges(links)

            # only if graph hasn't been seen before print it
            canonical = canonicalize_short(copy)
            if canonical in seen:
                continue

            seen.add(canonical)
            candidates.append("%d %d %s" % (kSys, kEnv, graph_to_adj(copy)))

        return candidates



    PATH_OUT = os.path.join(PATH_THISFILE, "out/")
    def print_candidates(graph, graph6_str, candidates):
        fn = os.path.join(PATH_OUT, filename_for_graph6(graph6_str) + ".adjm")
        with open(fn, "w") as f:
            f.writelines(c + "\n" for c in candidates)

        # mark graph as done and unlock
        filesystem_mark_done(graph6_str)
        filesystem_unlock(graph6_str)

    THREADS=2
    with multiprocessing.pool.Pool(processes=THREADS) as pool:
        for graph6_str, graph in generate_all_graphs(kSys):
            # did we already process this graph?
            if filesystem_is_done(graph6_str):
                continue

            # try locking
            if not filesystem_lock(graph6_str):
                continue

            # if successfully locked, run
            pool.apply_async(
                get_all_nonisomorphic_environments,
                (graph,),
                callback=partial(print_candidates, graph, graph6_str)
            )


        pool.close()
        pool.join()

    