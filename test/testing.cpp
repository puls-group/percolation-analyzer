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

TEST_CASE("The graph should detect one dimension per component", "[graph component]")
{

    size_t num_vertices = 30;

    TranslationVector trans0, trans1, trans2, trans3;

    trans0.vec[0] = 0;
    trans0.vec[1] = 0;
    trans0.vec[2] = 0;

    trans1.vec[0] = 1;
    trans1.vec[1] = 0;
    trans1.vec[2] = 0;

    trans2.vec[0] = 0;
    trans2.vec[1] = 1;
    trans2.vec[2] = 0;

    trans3.vec[0] = 0;
    trans3.vec[1] = 0;
    trans3.vec[2] = 1;

    SECTION("Not adding edges should yield one component per vertex")
    {
        PercolationGraph graph;
        graph.reserve_vertices(num_vertices);
        std::vector<ComponentInfo> components = graph.get_component_percolation_info();

        REQUIRE(components.size() == num_vertices);

        for (size_t i = 0; i < num_vertices; i++)
        {
            REQUIRE(components[i].component_index == i);
            REQUIRE(components[i].percolation_dim == 0);
            REQUIRE(components[i].vertices.size() == 1);
        }
    }

    PercolationGraph graph2;
    graph2.reserve_vertices(num_vertices);

    SECTION("Successively adding edges should reduce number of components")
    {

        for (size_t i = 0; i < num_vertices - 1; i++)
        {
            graph2.add_edge(i, i + 1, trans0);

            std::vector<ComponentInfo> components2 = graph2.get_component_percolation_info();

            REQUIRE(components2.size() == num_vertices - i - 1);
        }
    }

    SECTION("Creating a loop should still detect exactly one component")
    {

        for (size_t i = 0; i < num_vertices - 1; i++)
        {
            graph2.add_edge(i, i + 1, trans0);

            std::vector<ComponentInfo> components2 = graph2.get_component_percolation_info();

            REQUIRE(components2.size() == num_vertices - i - 1);
        }

        graph2.add_edge(num_vertices - 1, 0, trans0);

        std::vector<ComponentInfo> components3 = graph2.get_component_percolation_info();

        REQUIRE(components3.size() == 1);
    }
}

TEST_CASE("The graph should correctly detect 1 vertex loop dimensions", "[graph loop]")
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

    // 1D loopbacks
    graph.add_edge(0, 0, trans1);
    graph.add_edge(0, 0, trans1);
    graph.add_edge(1, 1, trans2);
    graph.add_edge(2, 2, trans3);

    // 2D loopbacks
    graph.add_edge(3, 3, trans1);
    graph.add_edge(3, 3, trans2);

    graph.add_edge(4, 4, trans1);
    graph.add_edge(4, 4, trans3);

    graph.add_edge(5, 5, trans2);
    graph.add_edge(5, 5, trans3);

    // 3D loopbacks
    graph.add_edge(6, 6, trans1);
    graph.add_edge(6, 6, trans2);
    graph.add_edge(6, 6, trans3);

    std::vector<ComponentInfo> components = graph.get_component_percolation_info();

    REQUIRE(components.size() == 7);

    // test 1d loops
    for (size_t i = 0; i < 3; i++)
    {
        REQUIRE(components[i].vertices.size() == 1);
        REQUIRE(components[i].percolation_dim == 1);
        REQUIRE(components[i].component_index == i);
    }
    // test 2d loops
    for (size_t i = 3; i < 6; i++)
    {
        REQUIRE(components[i].vertices.size() == 1);
        REQUIRE(components[i].percolation_dim == 2);
        REQUIRE(components[i].component_index == i);
    }
    // test 3d loops
    for (size_t i = 6; i < 7; i++)
    {
        REQUIRE(components[i].vertices.size() == 1);
        REQUIRE(components[i].percolation_dim == 3);
        REQUIRE(components[i].component_index == i);
    }
}

TEST_CASE("The graph should correctly detect 1 dimensional components", "[graph 1d]")
{
    PercolationGraph graph;

    TranslationVector trans1, trans2, trans3;

    trans1.vec[0] = GENERATE(-1, 1);
    trans1.vec[1] = GENERATE(0, 1);
    trans1.vec[2] = GENERATE(-1, 0);

    trans2 = trans1 + trans1;
    trans3 = -trans1;
    //std::cerr << trans1.vec[0] << "|" << trans1.vec[1] << "|" << trans1.vec[0] << std::endl;
    //std::cerr << trans2.vec[0] << "|" << trans2.vec[1] << "|" << trans2.vec[0] << std::endl;
    //std::cerr << trans3.vec[0] << "|" << trans3.vec[1] << "|" << trans3.vec[0] << std::endl;

    graph.add_edge(0, 1, trans1);
    graph.add_edge(1, 0, trans1);
    graph.add_edge(2, 3, trans2);
    graph.add_edge(3, 2, trans2);
    graph.add_edge(4, 5, trans3);
    graph.add_edge(5, 4, trans3);

    std::vector<ComponentInfo> components = graph.get_component_percolation_info();

    REQUIRE(components.size() == 3);

    if (trans1.vec[0] == 0 && trans1.vec[1] == 0 && trans1.vec[2] == 0)
    {
        for (size_t i = 0; i < 3; i++)
        {
            REQUIRE(components[i].vertices.size() == 2);
            REQUIRE(components[i].percolation_dim == 0);
            REQUIRE(components[i].component_index == i);
        }
    }
    else
    {
        for (size_t i = 0; i < 3; i++)
        {
            REQUIRE(components[i].vertices.size() == 2);
            REQUIRE(components[i].percolation_dim == 1);
            REQUIRE(components[i].component_index == i);
        }
    }
}

TEST_CASE("Deal with more complex 1d situations", "[graph 1d complex]")
{
    PercolationGraph graph;

    TranslationVector trans1, trans0, trans1_;

    trans1.vec[0] = 1;
    trans1.vec[1] = 0;
    trans1.vec[2] = 0;

    trans0.vec[0] = 0;
    trans0.vec[1] = 0;
    trans0.vec[2] = 0;

    trans1_ = -trans1;

    // 1D fan out
    graph.add_edge(0, 1, trans1);
    graph.add_edge(0, 2, trans1);
    graph.add_edge(0, 3, trans1);

    // Collect the fan
    graph.add_edge(1, 4, trans0);
    graph.add_edge(2, 4, trans0);
    graph.add_edge(3, 4, trans0);

    // Revert to another cell
    graph.add_edge(4, 5, trans1_);
    graph.add_edge(5, 6, trans1_);

    // Another fan out
    graph.add_edge(6, 7, trans0);
    graph.add_edge(6, 8, trans0);
    graph.add_edge(6, 9, trans0);
    graph.add_edge(6, 10, trans0);

    // Go to origin
    graph.add_edge(7, 0, trans1);
    // Will not do the 8
    graph.add_edge(9, 0, trans1);
    graph.add_edge(10, 0, trans1);

    std::vector<ComponentInfo> components = graph.get_component_percolation_info();

    REQUIRE(components.size() == 1);

    const ComponentInfo &component = components[0];

    REQUIRE(component.vertices.size() == 11);
    REQUIRE(component.percolation_dim == 0);
    REQUIRE(component.component_index == 0);
}