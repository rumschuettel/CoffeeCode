#!/usr/bin/env python3
# needs networkx, igraph
# as well as geng and vcolg from nauty27

## GRAPH
import networkx
import igraph
import base64

def to_graph6(graph, with_header=False):
    "igraph to graph6 ascii-encoded bytes"
    edges = graph.get_edgelist()
    nx_graph = networkx.Graph(edges)
    return networkx.readwrite.to_graph6_bytes(nx_graph, header=with_header).strip().decode("ascii")

def from_graph6(graph6_str):
    "graph6 ascii-encoded bytes to igraph"
    nx_graph = networkx.readwrite.from_graph6_bytes(b">>graph6<<" + graph6_str.strip().encode("ascii"))
    edges = list(nx_graph.edges)
    return igraph.Graph(edges)

def to_adjm(graph):
    "return adjacency matrix as a flattened string of 0s and 1s"
    return "".join([ str(i) for r in graph.get_adjacency() for i in r ])

def filename_for_graph6(graph6_str):
    return base64.b32encode(graph6_str.encode("ascii")).decode("ascii")



import itertools

def unordered_tuples(iterable, min_num_el, max_num_el):
    "tuples([1,2,3], 1, 2) --> (1,) (2,) (3,) (1,1) (1,2) (1,3) (2,2) (2,3) (3,3)"
    return itertools.chain.from_iterable(
        itertools.combinations_with_replacement(iterable, r) for r in range(min_num_el, max_num_el+1)
    )


def all_possible_vertex_subsets(kSys, kEnv):
    "iterator over all possible lists of system vertex subsets up to length kEnv"
    vsys = list(itertools.product(*[[0,1]]*kSys))[1:]
    return unordered_tuples(vsys, 1, kEnv)

def all_possible_environment_links(kSys, kEnv):
    "iterator over all possible system-environment edge lists"
    for subsets in all_possible_vertex_subsets(kSys, kEnv):
        yield ([
                (v_idx, subset_idx+kSys) # edges
                for subset_idx, subset in enumerate(subsets)
                for v_idx, v in enumerate(subset) if v==1
            ],
            len(list(subsets)) # kEnv
        )


## IO
import subprocess
import os

PATH_THISFILE = os.path.dirname(os.path.realpath(__file__))
PATH_GENG = os.path.join(PATH_THISFILE, "bin/geng")
PATH_SHORTG = os.path.join(PATH_THISFILE, "bin/shortg")

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

def run(cmd, *args):
    subprocess.run([cmd, *args])


## POOR MAN'S FILE SYSTEM LOCKS
PATH_FILESYSTEM_LOCK = os.path.join(PATH_THISFILE, "locks/")

def filesystem_temp_fn(graph6_str, extension):
    return os.path.join(PATH_FILESYSTEM_LOCK, filename_for_graph6(graph6_str) + "." + extension)

def filesystem_lock(graph6_str):
    "check whether we already processed graph by touching a file indicator"
    lockfile = filesystem_temp_fn(graph6_str, "lock")
    if os.path.isfile(lockfile):
        return False

    open(lockfile, "w").close()
    return True

def filesystem_unlock(graph6_str):
    "removes lock file indicator"
    lockfile = filesystem_temp_fn(graph6_str, "lock")
    try:
        os.remove(lockfile)
        return True
    except:
        return False

def filesystem_is_done(graph6_str):
    "checks if that graph was processed"
    donefile = filesystem_temp_fn(graph6_str, "done")
    return os.path.isfile(donefile)

def filesystem_mark_done(graph6_str):
    "marks a graph as done"
    donefile = filesystem_temp_fn(graph6_str, "done")
    open(donefile, "w").close()

## NAUTY LINKS
def generate_all_graphs(vertex_count):
    for line in pipe(PATH_GENG, "-c", "-q", vertex_count):
        yield line, from_graph6(line)

def filter_graphs(filename, kSys):
    run(PATH_SHORTG, "-f%s" % "a"*kSys, filename)


## PROGRAM
import multiprocessing.pool
from functools import partial


if __name__ == "__main__":
    kSys = 5
    max_kEnv = kSys
    PATH_OUT = os.path.join(PATH_THISFILE, "out/")

    def generate_all_nonisomorphic_environments(graph, graph6_str):
        print("doing", graph6_str)
        filename = filesystem_temp_fn(graph6_str, "redundant-graphs")

        count_redundant = 0
        with open(filename, "w", encoding="ascii") as f:
            # add environment links and store as a flat file of the graphs
            for links, kEnv in all_possible_environment_links(kSys, max_kEnv):
                # create copy of graph with environment vertices added
                copy = graph.copy()
                copy.add_vertices(kEnv)
                copy.add_edges(links)

                f.write("%s\n" % to_graph6(copy))
                count_redundant += 1

        print("filtering", graph6_str)
        # call nauty's shortg to prune graph list
        filter_graphs(filename, kSys)

        print("writing output for", graph6_str)
        # read list back in line by line and write final output
        count_unique = 0
        with open(filename, "r", encoding="ascii") as f_in:
            with open(os.path.join(PATH_OUT, filename_for_graph6(graph6_str) + ".adjm"), "w", encoding="ascii") as f_out:
                for with_env_graph6_str in f_in:
                    graph = from_graph6(with_env_graph6_str)
                    kEnv = graph.vcount() - kSys

                    f_out.write("%d %d %s %s" % (
                        kSys,
                        kEnv,
                        with_env_graph6_str,
                        to_adjm(graph)
                    ))

                    count_unique += 1


        # print some stats
        print("done %s: pruned from %d to %d unique" % (graph6_str, count_redundant, count_unique))

    def mark_done(graph6_str):
        # mark graph as done and unlock
        filesystem_mark_done(graph6_str)
        filesystem_unlock(graph6_str)

#    THREADS=1
#    with multiprocessing.pool.Pool(processes=THREADS) as pool:
    for graph6_str, graph in generate_all_graphs(kSys):
        # did we already process this graph?
        if filesystem_is_done(graph6_str):
            continue

        # try locking
        if not filesystem_lock(graph6_str):
            continue

        generate_all_nonisomorphic_environments(graph, graph6_str)
        mark_done(graph6_str)

        # if successfully locked, run
#        pool.apply_async(
#            generate_all_nonisomorphic_environments,
#            (graph, graph6_str),
#            callback=partial(mark_done, graph6_str)
#        )


#        pool.close()
#        pool.join()

    