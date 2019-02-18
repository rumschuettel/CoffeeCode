# CoffeeCode

## Requirements

* gcc 8 or newer/msvc 2017 or newer
* boost 1.69.0 source
* nauty 2.7rc1

The cxx compiler must be installed, boost and nauty extracted to some folders, respectively.
Note down these directories for the build instructions below.

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

| Option             | Value to enter, or comment.                            |
|--------------------|--------------------------------------------------------|
| `CMAKE_BUILD_TYPE` | `Release`                                              |
| `PATH_BOOST`       | Path to boost source directory. Must contain `boost/`. |
| `PATH_NAUTY`       | Path to nauty source directory. Must contain `nauty.c` |
| `K_ENV`            | Number of environment vertices.                        |
| `K_SYS`            | Number of system vertices.                             |
| `SYMMETRIC_SOLVER` | `ON` or `OFF`.                                         |

If `SYMMETRIC_SOLVER` is `OFF`, the nauty and boost directories are irrelevant.

Now press **[c]** again, and then **[g]** to generate the cmake scripts and exit ccmake.

Execute `cmake .` to create the build scripts, then compile with `make`.

## Advanced Configuration

Advanced config options (such as a path to a custom compiler, or optimization options) can be accessed from within ccmake after pressing **[t]**.

## Specifying Problem Instances

### Non-Symmetric Case

The full solver is run. `K_SYS` and `K_ENV` have to be set to select the most optimal internal representation of vectors, but the program can then be run on any instance with these parameters.

Run `./CoffeeCode`, then specify as input an adjacency matrix as a list of `0`s and `1`s, e.g. `010 101 010` for a path graph of three vertices (the 2-rep code). Whitespace is ignored. Confirm with **[enter]**.

### Symmetric Case

The symmetric solver is run. For optimal performance, the code has to be compiled _once per problem instance_; it suffices, however, to configure `K_SYS` and `K_ENV` with ccmake once, and then vary the instance file; the latter is specified in the file `cc-instance.h`, which contains a struct with four variables:

```
struct graphstate_instance {
	using sgs = SGSTransversal<
		SGSGenerator<1, Group<
		Permutation<2, 1, 0>
		>>
		>;
	constexpr static AdjacencyMatrixT<4> adjacency_matrix{ {{0, 1, 0, 0}, {1, 0, 1, 1}, {0, 1, 0, 0}, {0, 1, 0, 0}} };
};
```

`adjacency_matrix{ ... }` contains the adjacency matrix as a list of lists, where each list has to be terminated by curly braces. It is a `K_TOT x K_TOT` matrix, where `K_TOT=K_SYS+K_ENV`.

`sgs` is the strong generating set transversal given as a list of `SGSGenerator<...>`s. Each `SGSGenerator` specifies the stabilizer point, e.g. `SGSGenerator<7, ...>`, which is a number from `0...K_SYS-1`. Its second parameter is a `Group<...>`, specified as a list of `Permutation<...>`s.

Each permutation is given in list form and must contain each number from `0...K_SYS-1` exactly once (this is checked at compile time). For reasons of C++'s syntax, the lists here are each enclosed in angular brackets. Whitespace is ignored.

For ease of use, create the `cc-instance.h` file with Mathematica.