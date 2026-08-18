cmake_minimum_required(VERSION 3.18)
# Must be set before project(): CUDA compiler detection would otherwise
# initialize CMAKE_CUDA_ARCHITECTURES to the nvcc default (52) and the
# "if(NOT ...)" guard below project() would never fire.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 89)
endif()
project(zdem LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(ZDEM_PRECISION "float" CACHE STRING "floating point precision: float|double")
string(TOUPPER ${ZDEM_PRECISION} ZDEM_PRECISION_UP)
if(ZDEM_PRECISION_UP STREQUAL "DOUBLE")
    set(ZDEM_DEFINE ZDEM_DOUBLE)
endif()

enable_testing()

add_library(zdem_core
    src/geometry/stl_io.cpp
    src/geometry/mesh_build.cpp
    src/host/config_io.cpp
    src/host/vtk_io.cpp
    src/host/sim_build.cpp
    src/physics/pp_contact.cpp
    src/physics/wall_contact.cpp
    src/physics/integrate.cpp
)
target_include_directories(zdem_core PUBLIC src)

add_executable(zdem_cpu
    src/main.cpp
)
target_link_libraries(zdem_cpu PRIVATE zdem_core)

add_executable(zdem_gpu
    src/gpu/dem_gpu.cu
    src/gpu/dem_state.cu
)
target_link_libraries(zdem_gpu PRIVATE zdem_core)
target_compile_definitions(zdem_gpu PRIVATE ${ZDEM_DEFINE})
# Static cudart links the static CRT (LIBCMT) while MSVC host code uses /MD,
# which triggers LNK4098 and mixes two CRTs in one process. Keep one CRT:
# shared cudart matches the /MD default of zdem_core host code.
set_property(TARGET zdem_gpu PROPERTY CUDA_RUNTIME_LIBRARY Shared)

add_executable(test_math
    tests/test_math.cpp
)
target_link_libraries(test_math PRIVATE zdem_core)
add_test(NAME test_math COMMAND test_math)

add_executable(test_geometry
    tests/test_geometry.cpp
)
target_link_libraries(test_geometry PRIVATE zdem_core)
add_test(NAME test_geometry COMMAND test_geometry)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(
        NAME smoke_configs
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tests/smoke/run_smoke_tests.py
                --exe $<TARGET_FILE:zdem_cpu>
                --root ${CMAKE_SOURCE_DIR}
    )
    add_test(
        NAME wall_contact_regression
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tests/smoke/check_wall_contact_regression.py
                --exe $<TARGET_FILE:zdem_cpu>
                --root ${CMAKE_SOURCE_DIR}
    )
    add_test(
        NAME wall_multiloop_regression
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tests/smoke/check_wall_multiloop_regression.py
                --exe $<TARGET_FILE:zdem_cpu>
                --root ${CMAKE_SOURCE_DIR}
    )
    # GPU smoke: run zdem_gpu on the same static-contact scenario, assert
    # contacts >= 1 at step 0. SKIPs (exit 0) when no NVIDIA device exists.
    add_test(
        NAME gpu_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tests/smoke/check_gpu_smoke.py
                --exe $<TARGET_FILE:zdem_gpu>
                --root ${CMAKE_SOURCE_DIR}
    )
endif()
