# CoffeeCode

## Requirements

* gcc 8 or newer/msvc 2017 or newer
* boost 1.69.0 source
* nauty 2.7rc1

The cxx compiler must be installed, boost and nauty extracted to some folders, respectively.
Note down these directories for the build instructions below.

I recommend building both the symmetric and full solver in two separate directories (e.g. under `build/release_symm` and `build/release_full`); the paths can then be set in the Mathematica interface which then calls the symmetric or full solver accordingly.

Due to design decisions, the symmetric solver has to be built _once per problem instance_. This is such that static polymorphism can be exploited in full within the compiler to speed up the program's execution; this goes under the assumption that for large graphs, compiling the code once adds negligible time overhead.

The full solver has to be built _once per combination of `K_SYS` and `K_ENV`_, and can then be executed with just an adjacency matrix; the full solver also features a command line interface.

## Build Configuration

The following steps are done to create a *Release* version of CoffeeCode.
Similar steps can be done to build a *Debug* version.

Open a terminal and `cd` into this directory. Run the following commands:

```
    mkdir -p build/release
    cd build/release
    ccmake ../..
```

The list of available commands is shown at the bottom of the screen.
Press **[c]**, and set the following variables by first selecting the field, pressing **[enter]**, editing, then pressing **[enter]** again.

| Option                    | Value to enter, or comment.                            |
|---------------------------|--------------------------------------------------------|
| `CMAKE_BUILD_TYPE`        | `Release`                                              |
| `PATH_BOOST`              | Path to boost source directory. Must contain `boost/`. |
| `PATH_NAUTY`              | Path to nauty source directory. Must contain `nauty.c` |
| `PATH_CUSTOM_CC_INSTANCE` | Path to custom symmetric instance. Set in Mathematica. |
| `SYMMETRIC_SOLVER`        | `ON` or `OFF`.                                         |
| `REDUCE_LAMBDA`           | `ON` or `OFF`.                                         |

If `SYMMETRIC_SOLVER` is `OFF`, the nauty directory is irrelevant.
`REDUCE_LAMBDA` does polynomial simplification within C++, which is faster than doing this externally and reduces the output size.

Now press **[c]** again, and then **[g]** to generate the cmake scripts and exit ccmake.

Execute `cmake .` to create the build scripts.

### SYMMETRIC SOLVER

The symmetric solver is ready to be used from Mathematica; simply set the path there accordingly.

### FULL SOLVER

The full solver Mathematica interface assumes that a set of binaries have been precompiled (i.e. one for each combination of `K_SYS` and `K_ENV` for, say, up to 20 vertices - higher numbers are unlikely to make much sense and will take days to run due to the sheer number of iterations necessary).

For this, an automatic build script has been added to the build directory; simply run

```
    ./build_all.sh
```
and watch the magic happen. Once done, set the path in the MM interface accordingly.

## Advanced Configuration

Advanced config options (such as a path to a custom compiler, or optimization options) can be accessed from within ccmake after pressing **[t]**.

## Specifying Problem Instances

The easiest way to interface CoffeeCode is via the Mathematica script found under `scripts/`. Don't forget to set the paths up in the first few lines of the notebook.

### Non-Symmetric Case

Run `./CoffeeCode`, then specify as input an adjacency matrix as a list of `0`s and `1`s, e.g. `010 101 010` for a path graph of three vertices (the 2-rep code). Whitespace is ignored. Confirm with **[enter]**.

### Symmetric Case

The symmetric solver is run. For optimal performance, the code has to be compiled _once per problem instance_; due to the way nauty has to be linked we need to set `K_SYS` and `K_ENV` both as command line arguments for make (as `make K_SYS=3 K_ENV=1`) as well as *matching* within the file `cc-instance-custom.h`, which has been added to the build directory.
Note: this file is included from within `cc-instance.h`, which contains a struct with four variables:

```
struct graphstate_instance {
	using sgs = SGSTransversal<
		SGSGenerator<1, Group<
		Permutation<2, 1, 0>
		>>
		>;
	constexpr static size_t k_sys = 3, k_env = 1;
	constexpr static AdjacencyMatrixT<4> adjacency_matrix{ {{0, 1, 0, 0}, {1, 0, 1, 1}, {0, 1, 0, 0}, {0, 1, 0, 0}} };
};
```

`k_sys` and `k_env` have to match the command line parameters given.

`adjacency_matrix{ ... }` contains the adjacency matrix as a list of lists, where each list has to be terminated by curly braces. It is a `K_TOT x K_TOT` matrix, where `K_TOT=K_SYS+K_ENV`.

`sgs` is the strong generating set transversal given as a list of `SGSGenerator<...>`s. Each `SGSGenerator` specifies the stabilizer point, e.g. `SGSGenerator<7, ...>`, which is a number from `0...K_SYS-1`. Its second parameter is a `Group<...>`, specified as a list of `Permutation<...>`s.

Each permutation is given in list form and must contain each number from `0...K_SYS-1` exactly once (this is checked at compile time). For reasons of C++'s syntax, the lists here are each enclosed in angular brackets. Whitespace is ignored.

For ease of use, create the `cc-instance-custom.h` file with Mathematica, which does all the annoying redundant specifications automatically.