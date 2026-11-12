cmake_minimum_required(VERSION 3.16)
project(tacc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

add_library(graph_lib STATIC src/common/graph.cpp)
add_library(core_lib STATIC src/core/bicore_decomposition.cpp)
target_sources(core_lib PRIVATE src/core/dynamic_bicore.cpp)
add_library(index_lib STATIC
    src/index/tab_index.cpp
    src/index/vat_index.cpp
    src/index/uat_index.cpp
    src/index/mat_index.cpp
    src/index/index_builder.cpp
    src/index/index_io.cpp
)
add_library(query_lib STATIC
    src/query/online_query.cpp
    src/query/index_query_utils.cpp
    src/query/vat_query.cpp
    src/query/uat_query.cpp
    src/query/mat_query.cpp
    src/query/mat_recovery.cpp
)
add_library(experiments_lib STATIC
    src/experiments/experiment.cpp
    src/experiments/dataset_stats.cpp
    src/experiments/query_batch.cpp
    src/experiments/dataset_catalog.cpp
    src/experiments/workload.cpp
    src/experiments/online_baseline.cpp
    src/experiments/mat_maintenance.cpp
    src/experiments/experiment_cli.cpp
)

target_include_directories(graph_lib PUBLIC include)
target_include_directories(core_lib PUBLIC include)
target_include_directories(index_lib PUBLIC include)
target_include_directories(query_lib PUBLIC include)
target_include_directories(experiments_lib PUBLIC include)
target_link_libraries(core_lib PUBLIC graph_lib)
target_link_libraries(index_lib PUBLIC core_lib graph_lib)
target_link_libraries(query_lib PUBLIC index_lib core_lib graph_lib)
add_library(joint_index_lib STATIC
    src/index/joint_index.cpp
    src/index/joint_index_builder.cpp
    src/index/joint_index_opt.cpp
)
target_include_directories(joint_index_lib PUBLIC include)
target_link_libraries(joint_index_lib PUBLIC query_lib)
target_link_libraries(experiments_lib PUBLIC query_lib index_lib core_lib graph_lib)
target_link_libraries(experiments_lib PUBLIC joint_index_lib)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(experiments_lib PUBLIC stdc++fs)
    target_link_libraries(joint_index_lib PUBLIC stdc++fs)
endif()

add_executable(tacc main.cpp)
target_link_libraries(tacc PRIVATE experiments_lib)

option(BUILD_VAT_DIRECTORY_TOOL "Build the standalone VAT directory supplement tool" ON)
if(BUILD_VAT_DIRECTORY_TOOL)
    add_executable(tacc_vat_add_directory tools/vat_add_directory.cpp)
    target_link_libraries(tacc_vat_add_directory PRIVATE joint_index_lib)
endif()

include(CTest)
if(BUILD_TESTING)
    if(UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_executable(tacc_batch_query_timing_test tests/test_batch_query_timing.cpp)
        target_link_libraries(tacc_batch_query_timing_test PRIVATE joint_index_lib)
        target_link_options(tacc_batch_query_timing_test PRIVATE
            "-Wl,--wrap=_ZNSt6chrono3_V212steady_clock3nowEv")
        add_test(NAME batch_query_timing
            COMMAND tacc_batch_query_timing_test ${CMAKE_CURRENT_SOURCE_DIR}/data/small)
        set_tests_properties(batch_query_timing PROPERTIES TIMEOUT 120)
    endif()
    add_executable(tacc_mat_construction_test tests/test_mat_construction.cpp)
    target_link_libraries(tacc_mat_construction_test PRIVATE query_lib)
    add_test(NAME mat_construction_correctness
        COMMAND tacc_mat_construction_test ${CMAKE_CURRENT_SOURCE_DIR}/data/small)
    set_tests_properties(mat_construction_correctness PROPERTIES TIMEOUT 120)
    add_executable(tacc_joint_construction_test tests/test_joint_construction.cpp)
    target_link_libraries(tacc_joint_construction_test PRIVATE joint_index_lib)
    add_test(NAME joint_construction_correctness
        COMMAND tacc_joint_construction_test ${CMAKE_CURRENT_SOURCE_DIR}/data/small)
    set_tests_properties(joint_construction_correctness PROPERTIES TIMEOUT 120)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME exp8_dblp
            COMMAND ${Python3_EXECUTABLE} -B
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_exp8_dblp.py
                    $<TARGET_FILE:tacc>)
        set_tests_properties(exp8_dblp PROPERTIES TIMEOUT 120)
        add_test(NAME experiment_entrypoints
            COMMAND ${Python3_EXECUTABLE} -B
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_experiment_cli.py
                    $<TARGET_FILE:tacc>)
        set_tests_properties(experiment_entrypoints PROPERTIES TIMEOUT 240)
        add_test(NAME dataset_label_variants
            COMMAND ${Python3_EXECUTABLE} -B
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_prepare_real_datasets.py)
        set_tests_properties(dataset_label_variants PROPERTIES TIMEOUT 120)
        add_test(NAME mat_maintenance
            COMMAND ${Python3_EXECUTABLE} -B
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mat_maintenance.py
                    $<TARGET_FILE:tacc>)
        set_tests_properties(mat_maintenance PROPERTIES TIMEOUT 240)
        add_test(NAME exp5_joint_scalability
            COMMAND ${Python3_EXECUTABLE} -B
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_exp5_joint_scalability.py
                    $<TARGET_FILE:tacc>)
        set_tests_properties(exp5_joint_scalability PROPERTIES TIMEOUT 120)
        add_test(NAME balanced_query
            COMMAND ${Python3_EXECUTABLE} -B
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_balanced_query.py
                    $<TARGET_FILE:tacc>)
        set_tests_properties(balanced_query PROPERTIES TIMEOUT 120)
    endif()
endif()
