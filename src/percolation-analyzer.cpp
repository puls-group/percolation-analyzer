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

#include "percolation-analyzer.h"
#include "percolation-detection.hpp"

extern "C"
{
    struct PercolationAnalyzer
    {
        percolation::PercolationGraph graph;
    };

    struct PercolationAnalyzer *create_PercolationAnalyzer()
    {
        struct PercolationAnalyzer *res = new struct PercolationAnalyzer();
        return res;
    }

    void free_PercolationAnalyzer(struct PercolationAnalyzer *ptr)
    {
        if (ptr != nullptr)
            delete ptr;
    }

    void add_edge(struct PercolationAnalyzer *ptr, size_t source_node, size_t head_node, struct EdgeData *edge)
    {
        if (ptr == nullptr || edge == nullptr)
        {
            return;
        }
        percolation::PercolationGraph &graph = ptr->graph;
        percolation::EdgeData ed;
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            ed.translation[i] = edge->vec.pos[i];
        }

        graph.add_edge(source_node, head_node, ed);
    }

    void add_edge_translation(struct PercolationAnalyzer *ptr, size_t source_node, size_t head_node, struct TranslationVector *trans)
    {
        if (ptr == nullptr || trans == nullptr)
        {
            return;
        }
        percolation::PercolationGraph &graph = ptr->graph;
        percolation::TranslationVector tv;
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            tv[i] = trans->pos[i];
        }

        graph.add_edge(source_node, head_node, tv);
    }

    void reserve_memory(struct PercolationAnalyzer *ptr, size_t num_vertices)
    {
        if (ptr != nullptr)
            ptr->graph.reserve_vertices(num_vertices);
    }

    void set_vertex_data(struct PercolationAnalyzer *ptr, size_t vertex_index, struct VertexData *data)
    {
        if (ptr == nullptr || data == nullptr)
        {
            return;
        }
        percolation::PercolationGraph &graph = ptr->graph;
        percolation::VertexData vd;
        vd.index = data->index;

        graph.add_vertex(vertex_index, vd);
    }

    struct PercolationInfo
    {
        std::vector<percolation::ComponentInfo> components;
        size_t num_comps;
        size_t next_component;
    };

    struct PercolationInfo *getPercolationInfo(struct PercolationAnalyzer *ptr)
    {
        if (ptr == nullptr)
            return nullptr;

        struct PercolationInfo *res = new struct PercolationInfo;

        res->components = ptr->graph.get_component_percolation_info();
        res->num_comps = res->components.size();
        res->next_component = 0;
        return res;
    }

    void freePercolationInfo(struct PercolationInfo *ptr)
    {
        if (ptr != nullptr)
        {
            ptr->components.clear();
            delete ptr;
        }
    }

    size_t get_num_components(struct PercolationInfo *ptr)
    {
        if (ptr == nullptr)
            return 0;

        return ptr->num_comps;
    }

    struct ComponentInfo
    {
    };

    struct ComponentInfo *get_next_component(struct PercolationInfo *ptr)
    {
        if (ptr == nullptr)
            return nullptr;

        if (ptr->next_component >= ptr->num_comps)
        {
            return nullptr;
        }

        return (struct ComponentInfo *)&ptr->components[ptr->next_component++];
    }

    size_t get_component_index(struct ComponentInfo *ptr)
    {
        if (ptr == nullptr)
            return -1;

        return ((percolation::ComponentInfo *)ptr)->component_index;
    }

    size_t get_component_percolation_dimension(struct ComponentInfo *ptr)
    {
        if (ptr == nullptr)
            return 0;

        return ((percolation::ComponentInfo *)ptr)->percolation_dim;
    }

    struct VertexList
    {
        std::vector<percolation::VertexData> &verts;
        size_t num_vertices;
        size_t next_vertex;

        VertexList(std::vector<percolation::VertexData> &list) : verts(list) {}
    };

    struct VertexList *get_component_vertices(struct ComponentInfo *ptr)
    {
        if (ptr == nullptr)
            return nullptr;

        struct VertexList *res = new struct VertexList(((percolation::ComponentInfo *)ptr)->vertices);
        res->num_vertices = res->verts.size();
        res->next_vertex = 0;
        return res;
    }

    void free_component_vertices(struct VertexList *ptr)
    {
        if (ptr != nullptr)
        {
            delete ptr;
        }
    }

    int get_next_VertexData(struct VertexList *ptr, struct VertexData *data)
    {
        if (ptr == nullptr)
        {
            return 0;
        }

        if (ptr->next_vertex >= ptr->num_vertices)
        {
            return 0;
        }

        if (data != nullptr)
        {
            percolation::VertexData &vertdata = ptr->verts[ptr->next_vertex++];
            data->index = vertdata.index;
        }
        return 1;
    }

    size_t get_num_vertices(struct VertexList *ptr)
    {
        if (ptr == nullptr)
        {
            return 0;
        }
        return ptr->num_vertices;
    }
};