cmake_minimum_required(VERSION 3.14)
project(zdem_cpu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_library(zdem_core
    src/geometry/stl_io.cpp
    src/geometry/mesh_build.cpp
    src/host/config_io.cpp
    src/host/vtk_io.cpp
)
target_include_directories(zdem_core PUBLIC src)

add_executable(zdem_cpu
    src/main.cpp
)
target_link_libraries(zdem_cpu PRIVATE zdem_core)

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
endif()
