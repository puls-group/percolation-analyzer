#include "molecular_graph.hpp"

namespace mol
{

    MolecularGraph::MolecularGraph() : MolecularGraph(0) {}
    MolecularGraph::MolecularGraph(size_t num_atoms)
    {
        this->set_atom_count(num_atoms);
    }

    void MolecularGraph::set_atom_count(size_t num_atoms)
    {
        this->atom_positions.resize(num_atoms);
        this->bonds.resize(num_atoms);
    }

    bool MolecularGraph::set_basis(const std::vector<vec<graph_precision_type>> &triclinic_basis)
    {
        if (triclinic_basis[0][1] != 0.0 || triclinic_basis[0][2] != 0.0 || triclinic_basis[1][2] != 0.0 || triclinic_basis[0][0] == 0.0 || triclinic_basis[1][1] == 0.0 || triclinic_basis[2][2] == 0.0)
        {
            // Basis does not have triclinic properties
            return false;
        }
        this->triclinic_basis = triclinic_basis;
        return true;
    }

    bool MolecularGraph::set_atom_position(size_t atom_index, const vec<graph_precision_type> &pos)
    {
        if (atom_index >= n_atoms)
        {
            return false;
        }
        this->atom_positions[atom_index] = pos;
        return true;
    }

    bool MolecularGraph::add_bond(size_t atom_index_1, size_t atom_index_2)
    {

        if (atom_index_1 >= n_atoms || atom_index_2 >= n_atoms)
        {
            return false;
        }
        bonds[atom_index_1].push_back(atom_index_2);
        bonds[atom_index_2].push_back(atom_index_1);
        return true;
    }

    percolation::PercolationGraph MolecularGraph::get_percolation_graph() const
    {
        percolation::PercolationGraph res;

        return res;
    }
}