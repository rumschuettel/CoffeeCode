cmake_minimum_required (VERSION 3.13)

project (CoffeeCode)

message(STATUS "Info: K_SYS and K_ENV have to be passed via command line to make, e.g. 'K_SYS=3 K_ENV=1 make'.")
message(STATUS "      If not done, make will fail to build with a bunch of errors.")

option(SYMMETRIC_SOLVER "Compile symmetric solver" ON)

if (SYMMETRIC_SOLVER)
    set(PATH_CUSTOM_CC_INSTANCE "/home/jkrb2/programming/CoffeeCode/build/release/cc-instance-custom.h" CACHE PATH "File for custom symmetric CC instance")
    set(PATH_NAUTY "/home/jkrb2/programming/libs/nauty27rc1/" CACHE PATH "Nauty 2.7 RC1 path")
    set(PATH_BOOST "/home/jkrb2/programming/libs/boost/boost_1_69_0/" CACHE PATH "Boost 1.69 path")
endif()

set(CMAKE_CXX_COMPILER "/home/jkrb2/opt/anaconda5/envs/gcc8/bin/x86_64-conda_cos6-linux-gnu-g++" CACHE FILE "cxx compiler")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
add_compile_options(-march=native -mtune=generic -fmax-errors=5)
 
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
    )
    target_include_directories(Nauty
        PRIVATE
            ${PATH_NAUTY}
    )
    target_compile_definitions(Nauty
        PRIVATE
            "MAXN=(\$\(K_SYS\)+\$\(K_ENV\))"
    )
endif()


##### CoffeeCode
################
add_executable(CoffeeCode "CoffeeCode.cpp")
# $$K_SYS and $$K_ENV are processed from make from the environment,
# and they HAVE to be set!
target_compile_definitions(CoffeeCode PRIVATE K_SYS=\$\(K_SYS\) K_ENV=\$\(K_ENV\))
target_compile_options(CoffeeCode
    PRIVATE
        -Wall
        -Wextra
        -pedantic
        -pedantic-errors
        -Wconversion
        -Wsign-conversion
)

if (SYMMETRIC_SOLVER)
    target_include_directories(CoffeeCode
        SYSTEM PRIVATE
            ${PATH_NAUTY}
            ${PATH_BOOST}
        PRIVATE
            ${PROJECT_BINARY_DIR}
    )
    target_compile_definitions(CoffeeCode
        PRIVATE
            "MAXN=(K_SYS+K_ENV)"
            SYMMETRIC_SOLVER
    )
    target_sources(CoffeeCode
        PRIVATE
            CCSymmetric.cpp
    )
    target_link_libraries(CoffeeCode Nauty)

    configure_file(
        "${PROJECT_SOURCE_DIR}/cc-instance.h.in"
        "${PROJECT_BINARY_DIR}/cc-instance.h"
    )
else()
    target_sources(CoffeeCode PRIVATE "CCFull.cpp")
endif()