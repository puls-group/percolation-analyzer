#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "percolation-detection.hpp"

using namespace percolation;

TEST_CASE("The graph should automatically allocate memory for entries", "[graph allocation]")
{
    PercolationGraph graph;

    TranslationVector trans1, trans2, trans3;

    trans1.vec[0] = 1;
    trans1.vec[1] = 0;
    trans1.vec[2] = 0;

    trans2.vec[0] = 0;
    trans2.vec[1] = 1;
    trans2.vec[2] = 0;

    trans3.vec[0] = 0;
    trans3.vec[1] = 0;
    trans3.vec[2] = 1;

    graph.add_edge(0, 1, trans1);
    graph.add_edge(0, 2, trans1);
    graph.add_edge(0, 2, trans1);
}