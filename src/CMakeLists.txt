
add_library(percolation-analyzer-cpp percolation-detection.cpp)
target_include_directories(percolation-analyzer-cpp PUBLIC ${INCLUDE_DIR})

add_library(molecular-graph-cpp molecular_graph.cpp)
target_include_directories(molecular-graph-cpp PUBLIC ${INCLUDE_DIR})