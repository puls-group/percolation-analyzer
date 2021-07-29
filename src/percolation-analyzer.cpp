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