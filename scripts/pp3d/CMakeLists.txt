cmake_minimum_required( VERSION 3.0 )
project( pp3d )

add_compile_options(-march=native -mtune=native -fmax-errors=5 -fext-numeric-literals)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable( pp3d ./pp3d.cpp )

find_package(nlohmann_json 3.2.0 REQUIRED)
target_link_libraries(pp3d PRIVATE nlohmann_json::nlohmann_json)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories( ~/opt/anaconda5/envs/gcc8/include )
target_link_libraries(pp3d PRIVATE quadmath)