cmake_minimum_required (VERSION 3.5)

project (CoffeeCode)

set(K_SYS 3 CACHE STRING "System vertices")
set(K_ENV 1 CACHE STRING "Environment vertices")
option (SYMMETRIC_SOLVER "Compile symmetric solver" ON)
if (SYMMETRIC_SOLVER)
    set(PATH_NAUTY "/home/jkrb2/programming/libs/nauty27rc1/" CACHE PATH "Nauty 2.7 RC1 path")
    set(PATH_BOOST "/home/jkrb2/programming/libs/boost/boost_1_69_0/" CACHE PATH "Boost 1.69 path")
endif()

set(CMAKE_CXX_COMPILER "/home/jkrb2/opt/anaconda5/envs/gcc8/bin/x86_64-conda_cos6-linux-gnu-g++" CACHE FILE "cxx compiler")

 
##### CoffeeCode
add_executable(CoffeeCode "CoffeeCode.cpp")

target_compile_definitions(CoffeeCode PRIVATE K_SYS=${K_SYS} K_ENV=${K_ENV})

# add source files
if (SYMMETRIC_SOLVER)
    target_include_directories(CoffeeCode PRIVATE "${PATH_NAUTY}" "${PATH_BOOST}")
    target_compile_definitions(CoffeeCode PRIVATE "MAXN=(K_SYS+K_ENV)" SYMMETRIC_SOLVER)

    target_sources(CoffeeCode PRIVATE "CCSymmetric.cpp")
    target_sources(CoffeeCode PRIVATE "${PATH_NAUTY}/nauty.c" "${PATH_NAUTY}/nautil.c" "${PATH_NAUTY}/naugraph.c" "${PATH_NAUTY}/schreier.c" "${PATH_NAUTY}/naurng.c")

    configure_file("${PROJECT_SOURCE_DIR}/cc-instance.h" "${PROJECT_BINARY_DIR}/cc-instance.h" COPYONLY)
else()
    target_sources(CoffeeCode PRIVATE "CCFull.cpp")
endif()
