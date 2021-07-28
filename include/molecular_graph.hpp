#ifndef __MOLECULAR_GRAPH_H__
#define __MOLECULAR_GRAPH_H__

#include <vector>
#include "vec.hpp"

#include "percolation-detection.hpp"

namespace mol
{
    using graph_precision_type = double;

    class MolecularGraph
    {
    public:
        MolecularGraph();
        MolecularGraph(size_t num_atoms);

        void set_atom_count(size_t num_atoms);
        bool set_basis(const std::vector<vec<graph_precision_type>> &triclinic_basis);

        bool set_atom_position(size_t atom_index, const vec<graph_precision_type> &pos);
        bool add_bond(size_t atom_index_1, size_t atom_index_2);

        percolation::PercolationGraph get_percolation_graph() const;

    protected:
        size_t n_atoms;
        std::vector<vec<graph_precision_type>> triclinic_basis;

        std::vector<vec<graph_precision_type>> atom_positions;
        std::vector<std::vector<size_t>> bonds;
    };
}

#endif