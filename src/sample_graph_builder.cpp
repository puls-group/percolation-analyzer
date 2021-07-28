#include <iostream>
#include <cstdlib>
#include <cmath>
#include <random>

#include "molecular_graph.hpp"

/**
 * @brief Sample program to randomly generate a molecular graph and convert it into the percolation graph
 * 
 * In this program, we will generate a random set of points to populate a cuboid space and generate random bonds between those to build a sample graph to run the percolation analyzer on.
 * This can be replaced by a reading in routine for your respective trajectory format to populate the graph as well as the links/bonds between the atoms by one frame's data.
 * It is mainly meant to provide a sample for how to convert the atomic position and bond data into the percolation graph
 * 
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, char *argv[])
{
    // These are parameters for the generation of the random graph.
    // If you read in your graph, you do not need these
    size_t num_points = 2000;

    double dimx = 10;
    double dimy = 10;
    double dimz = 10;

    double volume = dimx * dimy * dimz;

    double volume_per_node = volume / double(num_points);
    double max_cuttoff = 2.0 * std::pow(volume_per_node, 1. / 3.0);
    double link_probability_below_cutoff = 0.7;

    // Let us set up the c++ way to generate random numbers.
    std::random_device rd;

    //
    // We choose a Mersenne twister engine
    //
    std::mt19937 engine(rd());

    //
    // Distributions for random numbers
    //
    std::uniform_real_distribution<double> distx(0., dimx);
    std::uniform_real_distribution<double> disty(0., dimy);
    std::uniform_real_distribution<double> distz(0., dimz);

    std::uniform_real_distribution<double> link_distr(0., 1.0);

    //
    // If you read in your molecular graph from a trajectory, this is where to start setting up
    // The respective data structure:
    //
    mol::MolecularGraph mol_graph;

    // Reserve memory for the number of atoms in the frame
    mol_graph.set_atom_count(num_points);

    // Set up a cuboid basis set for the system
    // Can be made triclinic as long as basis[0][1] == basis[0][2] == basis[1][2] == 0
    // as well as basis[0][0] != 0,  basis[1][1] != 0, basis[2][2] != 0
    // The precise setup of this will depend on your chosen basis

    std::cout << "Register the system pbc" << std::endl;
    std::vector<vec<double>> basis(3);
    basis[0].x = dimx;
    basis[1].y = dimy;
    basis[2].z = dimz;

    // Check the return value of set_base to make sure that the graph basis is correctly formatted
    if (!mol_graph.set_basis(basis))
    {
        std::cerr << "Error setting the molecular graph basis. Basis is not triclinic." << std::endl;
        exit(1);
    }

    // Generate (or read in) the position information of the atoms
    std::vector<vec<double>> positions(num_points);

    const double sq_cutoff = max_cuttoff * max_cuttoff;
    std::cout << "Randomly populate molecular graph" << std::endl;
    //
    // If you read in your graph data, this is where you populate the graph with your data instead of randomly generated positions and links.
    //
    for (size_t i = 0; i < num_points; i++)
    {
        positions[i].x = distx(engine);
        positions[i].y = disty(engine);
        positions[i].z = distz(engine);

        mol_graph.set_atom_position(i, positions[i]);

        for (size_t prev = 0; prev < i; prev++)
        {
            vec<double> diff = positions[prev] - positions[i];

            // Deal with the periodicity:
            vec<double> basis_coeffficients = decompose(diff, basis);
            // Choose the shortest connection as the basis by limiting the coefficients to [-0.5,0.5]

            // Calculate the shortest connecting vector respecting pbc and its square length
            vec<double> per_diff = basis[0] * basis_coeffficients[0] + basis[1] * basis_coeffficients[1] + basis[2] * basis_coeffficients[2];
            double distance2 = per_diff.norm2();

            // Randomly add bonds if below cutoff
            if (distance2 < sq_cutoff && link_distr(engine) <= link_probability_below_cutoff)
            {
                mol_graph.add_bond(prev, i);
            }
        }
    }

    //
    // At this point, mol_graph holds the entire graph data.
    // We will now convert it to the percolation detection graph
    //

    std::cout << "Convert molecular graph to percolation graph" << std::endl;
    // Retrieve the Percolation Graph associated with the molecular graph
    percolation::PercolationGraph percolation_graph = mol_graph.get_percolation_graph();

    std::cout << "Retrieve percolation data" << std::endl;
    // Now calculate the percolation info
    std::vector<percolation::ComponentInfo> percolation_data = percolation_graph.get_component_percolation_info();

    // Let us output the percolation info

    std::cout << "Results of the percolation analysis:" << std::endl;

    std::cout << "The system has a total of " << percolation_data.size() << " disjoint molecules" << std::endl;
    std::cout << "We will now list the data for each of these molecules..." << std::endl
              << std::endl;

    for (size_t c_index = 0; c_index < percolation_data.size(); c_index++)
    {
        std::cout << "Molecule #" << c_index << " stats:" << std::endl;

        std::cout << "Percolation dimension: " << percolation_data[c_index].percolation_dim << std::endl;
        std::cout << "#Atoms linked: " << percolation_data[c_index].vertices.size() << std::endl;

        // If you want full info on which atoms are linked, uncomment the following:
        /*auto &c_verts = percolation_data[c_index].vertices;
        std::cout << "Indices of atoms in this molecule:" << std::endl;
        for (size_t v_index = 0; v_index < c_verts.size(); v_index++)
        {
            std::cout << "\t" << c_verts[v_index].index << std::endl;
        }*/

        std::cout << std::endl;
    }
}