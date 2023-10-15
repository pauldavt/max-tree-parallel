cmake_minimum_required(VERSION 2.8.12)

project(pmt)

add_library(pmt STATIC
  src/misc/bit_array.cc
  src/misc/logger.cc
  src/parallel/thread_pool.cc
)

target_include_directories(pmt PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

set(CMAKE_CXX_FLAGS  "-Wfatal-errors -O3 -march=native -pedantic -std=c++14 ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_EXTENSIONS OFF)

link_libraries(pmt)

add_executable(trie_queue tests/trie_queue.cc)
add_executable(bit_array tests/bit_array.cc)

link_libraries(pthread atomic)

add_executable(parallel_for tests/parallel_for.cc)
add_executable(radix_sort tests/radix_sort.cc)
add_executable(image_blocks_2d tests/image_blocks_2d.cc)
add_executable(image_blocks_3d tests/image_blocks_3d.cc)
add_executable(maxtree_2d tests/maxtree_2d.cc)
add_executable(maxtree_3d tests/maxtree_3d.cc)
add_executable(connected_components tests/connected_components.cc)
add_executable(euler_tour_scan tests/euler_tour_scan.cc)
add_executable(tree_scan tests/tree_scan.cc)
add_executable(rootfix tests/rootfix.cc)
add_executable(direct_filter tests/direct_filter.cc)

add_executable(area_opening area_opening.cc)

target_include_directories(area_opening PUBLIC
  /usr/include/opencv4
)

target_link_libraries(area_opening opencv_core opencv_imgcodecs)
