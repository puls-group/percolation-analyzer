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

#ifndef __PERCOLATION_ANALYZER_H__
#define __PERCOLATION_ANALYZER_H__

#include <stddef.h>

extern "C"
{
    struct PercolationAnalyzer;

    struct TranslationVector
    {
        union
        {
            struct
            {
                long long x, y, z;
            };
            long long pos[3];
        };
    };

    struct EdgeData
    {
        struct TranslationVector vec;
    };

    struct VertexData
    {
        size_t index;
    };

    struct PercolationInfo;
    struct ComponentInfo;
    struct VertexList;

    struct PercolationAnalyzer *create_PercolationAnalyzer();

    void free_PercolationAnalyzer(struct PercolationAnalyzer *);

    void add_edge(struct PercolationAnalyzer *, size_t source_node, size_t head_node, struct EdgeData *);
    void add_edge_translation(struct PercolationAnalyzer *, size_t source_node, size_t head_node, struct TranslationVector *);

    void reserve_memory(struct PercolationAnalyzer *, size_t num_vertices);

    void set_vertex_data(struct PercolationAnalyzer *, size_t vertex_index, struct VertexData *);

    struct PercolationInfo *get_PercolationInfo(struct PercolationAnalyzer *);
    void free_PercolationInfo(struct PercolationInfo *);

    size_t get_num_components(struct PercolationInfo *);
    struct ComponentInfo *get_next_component(struct PercolationInfo *);

    size_t get_component_index(struct ComponentInfo *);
    size_t get_component_percolation_dimension(struct ComponentInfo *);
    struct VertexList *get_component_vertices(struct ComponentInfo *);

    void free_component_vertices(struct VertexList *);

    int get_next_VertexData(struct VertexList *, struct VertexData *);
    size_t get_num_vertices(struct VertexList *);
};

#endif