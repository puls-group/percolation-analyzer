
#include "percolation-detection.hpp"
#include <set>
namespace percolation
{

    struct TranslationVector TranslationVector::operator-() const
    {
        TranslationVector res;
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            res.vec[i] = -vec[i];
        }
    }

    struct TranslationVector TranslationVector::operator+(struct TranslationVector const &other) const
    {
        TranslationVector res;
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            res.vec[i] = vec[i] + other.vec[i];
        }
    }

    struct EdgeData EdgeData::inverse() const
    {
        EdgeData res(*this);
        res.translation = -this->translation;
        return res;
    }

    bool check_translation_independent(const std::vector<TranslationVector> &existing_base, const TranslationVector &new_vector)
    {
        return false;
    }

    bool PercolationGraph::reserve_vertices(size_t max_index)
    {
        if (this->vertices.size() > max_index)
        {
            return true;
        }
        this->vertices.resize(max_index);
        this->edges.resize(max_index);
        return true;
    }

    bool PercolationGraph::add_vertex(size_t vertex_index, const VertexData &vertex_data)
    {
        if (!reserve_vertices(vertex_index))
        {
            return false;
        }
        vertices[vertex_index] = vertex_data;
        return true;
    }

    bool PercolationGraph::add_edge(size_t vertex_index_base, size_t vertex_index_head, const EdgeData &edge_data)
    {
        size_t max_index = (vertex_index_base > vertex_index_head ? vertex_index_base : vertex_index_head);
        if (!reserve_vertices(max_index))
        {
            return false;
        }
        edges[vertex_index_base].push_back({vertex_index_head, edge_data});
        edges[vertex_index_head].push_back({vertex_index_base, edge_data.inverse()});
        return true;
    }

}