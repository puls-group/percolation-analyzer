/*
 * SPDX-FileCopyrightText: 2020 Kevin Höllring for PULS Group <kevin.hoellring@fau.de>
 *
 * SPDX-License-Identifier: MIT
 * 
 * Copyright (c) 2020 Kevin Höllring for PULS Group <kevin.hoellring@fau.de>
 * 
 * Authors: 2020 Kevin Höllring <kevin.hoellring@fau.de>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the “Software”), to deal 
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
 * copies of the Software, and to permit persons to whom the Software is 
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice (including the next paragraph) 
 * shall be included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE. 
 */

#include "molecular-graph.hpp"

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
        n_atoms = num_atoms;
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
        // Resize the percolation graph appropriately
        res.reserve_vertices(n_atoms);

        // Transform all positions into normalized basis components
        // This is equivalent to moving them all into one pbc cell and transforming the coordinates
        // into cuboid shape, which makes everything simpler.
        std::vector<vec<graph_precision_type>> normalized_positions(n_atoms);
        for (size_t base = 0; base < n_atoms; base++)
        {
            normalized_positions[base] = normalize_basis_coefficients(decompose(atom_positions[base], triclinic_basis));
        }

        // Parse the edges to be added:
        for (size_t base = 0; base < n_atoms; base++)
        {
            const vec<graph_precision_type> &b_pos = normalized_positions[base];
            for (size_t n = 0; n < bonds[base].size(); n++)
            {
                size_t head = bonds[base][n];
                if (head < base)
                {
                    // Both directions of the edge were created, we only need to consider one.
                    continue;
                }

                const vec<graph_precision_type> &h_pos = normalized_positions[head];

                // Let us build the correct translation vector
                percolation::TranslationVector trans;

                for (size_t dim = 0; dim < vector_space_dimension; dim++)
                {
                    if (b_pos[dim] < h_pos[dim])
                    {
                        // The base lies below the head in the cell
                        if (h_pos[dim] - b_pos[dim] > 0.5)
                        {
                            // The shorter connection intersects the cell boundary going down from the base
                            trans[dim] = -1;
                        }
                        else
                        {
                            // The connection is within the cell
                            trans[dim] = 0;
                        }
                    }
                    else
                    {
                        // Need to check edge crossing the other way around
                        // The base lies below the head in the cell
                        if (b_pos[dim] - h_pos[dim] > 0.5)
                        {
                            // The shorter connection intersects the cell boundary going up from the base
                            trans[dim] = +1;
                        }
                        else
                        {
                            // The connection is within the cell
                            trans[dim] = 0;
                        }
                    }
                }
                res.add_edge(base, head, trans);
            }
        }

        return res;
    }
}