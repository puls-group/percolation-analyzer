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