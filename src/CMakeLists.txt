# CPP interface for library
add_library(percolation-analyzer-cpp percolation-detection.cpp)
target_include_directories(percolation-analyzer-cpp PUBLIC ${INCLUDE_DIR})

# C wrapper for library
add_library(percolation-analyzer-c percolation-analyzer.cpp)
target_include_directories(percolation-analyzer-c PUBLIC ${INCLUDE_DIR})
target_link_libraries(percolation-analyzer-c percolation-analyzer-cpp)

# Cpp interface for molecular graph structure
add_library(molecular-graph-cpp molecular_graph.cpp)
target_include_directories(molecular-graph-cpp PUBLIC ${INCLUDE_DIR})