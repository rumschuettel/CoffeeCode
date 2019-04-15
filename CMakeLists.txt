cmake_minimum_required (VERSION 3.13)

project (CoffeeCode)

message(STATUS "Info: K_SYS and K_ENV have to be passed via command line to make, e.g. 'make K_SYS=3 K_ENV=1'.")
message(STATUS "      If not done, make will fail to build with a bunch of errors.")
message(STATUS "      Building the symmetric solver? MM will compile for you.")
message(STATUS "      Building the full solver? Run ./build-all.sh before using the MM interface.")

option(SYMMETRIC_SOLVER "Compile symmetric solver" ON)
option(REDUCE_LAMBDA "Merge lambda and lambda_a before printing" ON)
option(OPTIMIZE_FOR_DEPOLARIZING "Optimize polynomial computations for depolarizing channel" OFF)
option(FLOATING_POINT_MULTIPLICITY "Use approximate large numbers instead of exact integers" OFF)
option(PARALLELIZE "Use OpenMP to parallelize; number of threads picked automatically" OFF)

if (SYMMETRIC_SOLVER)
    set(PATH_CUSTOM_CC_INSTANCE "${PROJECT_BINARY_DIR}/cc-instance-custom.h" CACHE PATH "File for custom symmetric CC instance")
    set(PATH_NAUTY "/home/jkrb2/programming/libs/nauty27rc2/" CACHE PATH "Nauty 2.7 RC2 path")
    set(PATH_BOOST "/home/jkrb2/programming/libs/boost/boost_1_69_0/" CACHE PATH "Boost 1.69 path")
endif()

# empty odd defaults
set(CMAKE_CXX_FLAGS "")
set(CMAKE_C_FLAGS "")

add_compile_options(-march=native -fmax-errors=5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS OFF)


if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_compile_options(-diag-disable=all -prof-use)
endif()


##### Build Targets
###################

if (SYMMETRIC_SOLVER)
    ##### Nauty as separate library to link against, speeds up compilation
    ###########

    add_library(Nauty
        "${PATH_NAUTY}/nauty.c"
        "${PATH_NAUTY}/nautil.c"
        "${PATH_NAUTY}/naugraph.c"
        "${PATH_NAUTY}/schreier.c"
        "${PATH_NAUTY}/naurng.c"
        "${PATH_NAUTY}/naugroup.c"
        "${PATH_NAUTY}/gtnauty.c"
        "${PATH_NAUTY}/nausparse.c"
    )
    target_include_directories(Nauty
        PRIVATE
            ${PATH_NAUTY}
    )
    target_compile_definitions(Nauty
        PRIVATE
            "MAXN=(\$\(K_SYS\)+\$\(K_ENV\))"
            WORDSIZE=64
    )
endif()


##### CoffeeCode
################
add_executable(CoffeeCode "CoffeeCode.cpp")
# $$K_SYS and $$K_ENV are processed from make from the environment,
# and they HAVE to be set!
target_compile_definitions(CoffeeCode
    PRIVATE
        K_SYS=\$\(K_SYS\)
        K_ENV=\$\(K_ENV\)
)
if (REDUCE_LAMBDA)
    target_compile_definitions(CoffeeCode PRIVATE REDUCE_LAMBDA)
endif()
if (OPTIMIZE_FOR_DEPOLARIZING)
    target_compile_definitions(CoffeeCode PRIVATE OPTIMIZE_FOR_DEPOLARIZING)
endif()
if (FLOATING_POINT_MULTIPLICITY)
    target_compile_definitions(CoffeeCode PRIVATE FLOATING_POINT_MULTIPLICITY)
endif()
if (PARALLELIZE)
    target_compile_definitions(CoffeeCode PRIVATE PARALLELIZE)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(CoffeeCode PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

target_compile_options(CoffeeCode
    PRIVATE
        -Wall
        -Wextra
        -pedantic
        -pedantic-errors
        $<$<CONFIG:DEBUG>:-Wconversion>
        $<$<CONFIG:DEBUG>:-Wsign-conversion>
)
target_include_directories(CoffeeCode
    SYSTEM PRIVATE
        ${PATH_BOOST}
)

if (SYMMETRIC_SOLVER)
    target_include_directories(CoffeeCode
        SYSTEM PRIVATE
            ${PATH_NAUTY}
        PRIVATE
            ${PROJECT_BINARY_DIR}
    )
    target_compile_definitions(CoffeeCode
        PRIVATE
            "MAXN=(K_SYS+K_ENV)"
            WORDSIZE=64
            SYMMETRIC_SOLVER
    )
    target_sources(CoffeeCode
        PRIVATE
            CCSymmetric.cpp
            nauty-color-generator-vcolg.cpp
    )
    target_link_libraries(CoffeeCode PRIVATE Nauty)

    configure_file(
        "${PROJECT_SOURCE_DIR}/cc-instance.h.in"
        "${PROJECT_BINARY_DIR}/cc-instance.h"
    )
else()
    target_sources(CoffeeCode PRIVATE "CCFull.cpp")

    configure_file(
        "${PROJECT_SOURCE_DIR}/build-all.sh.in"
        "${PROJECT_BINARY_DIR}/build-all.sh"
    )
endif()
