#ifndef __MOLECULAR_GRAPH_H__
#define __MOLECULAR_GRAPH_H__

#include <vector>
#include "vec.hpp"

#include "percolation-detection.hpp"

namespace mol
{

    template <typename T>
    class MolecularGraph
    {
    public:
        MolecularGraph();
        MolecularGraph(size_t num_atoms);

        void set_atom_count(size_t num_atoms);
        bool set_basis(const std::vector<vec<T>> &triclinic_basis);

        void set_atom_position(size_t atom_index, const vec<T> &pos);
        void add_bond(size_t atom_index_1, size_t atom_index_2);

        percolation::PercolationGraph get_percolation_graph() const;

    protected:
        size_t n_atoms;
        std::vector<vec<T>> triclinic_basis;

        std::vector<vec<T>> atom_positions;
        std::vector<std::vector<size_t>> bonds;
    };
}

#endif