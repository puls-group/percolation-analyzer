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