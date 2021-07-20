
#include "percolation-detection.hpp"
#include <queue>
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
        return res;
    }

    struct TranslationVector TranslationVector::operator+(struct TranslationVector const &other) const
    {
        TranslationVector res;
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            res.vec[i] = vec[i] + other.vec[i];
        }
        return res;
    }

    struct TranslationVector TranslationVector::operator-(struct TranslationVector const &other) const
    {
        TranslationVector res;
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            res.vec[i] = vec[i] - other.vec[i];
        }
        return res;
    }

    bool TranslationVector::operator==(struct TranslationVector const &other) const
    {
        for(size_t i=0; i < vector_space_dimension; i++){
            if (vec[i] != other.vec[i]){
                return false;
            }
        }
        return true;
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
        vertices[vertex_index].index = vertex_index;

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


    std::vector<ComponentInfo> PercolationGraph::get_components() const
    {
        std::vector<ComponentInfo> component_info;
        size_t curr_comp = 0;
        size_t num_vertices = this->vertices.size();
        std::vector<bool> visited(num_vertices, false);

        for (size_t curr_vertex = 0; curr_vertex < num_vertices; curr_vertex++)
        {
            if (visited[curr_vertex])
            {
                continue;
            }

            std::queue<size_t, std::set<size_t>> vertex_queue;
            vertex_queue.push(curr_vertex);

            // Create new component metadata
            ComponentInfo new_comp;
            new_comp.component_index = curr_comp++;
            new_comp.perolation_dim = (size_t)-1;

            // Find all vertices in component via BFS
            while (!vertex_queue.empty())
            {
                size_t vert_index = vertex_queue.front();
                vertex_queue.pop();

                if (visited[vert_index])
                {
                    continue;
                }

                new_comp.vertices.push_back(vertices[vert_index]);
                visited[vert_index] = true;
                for (auto edge : this->edges[vert_index])
                {
                    size_t neighbor = edge.first;
                    if (!visited[neighbor])
                    {
                        vertex_queue.push(neighbor);
                    }
                }
            }
            component_info.push_back(new_comp);
        }
        return component_info;
    }
}