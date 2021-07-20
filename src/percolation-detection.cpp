
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
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            if (vec[i] != other.vec[i])
            {
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

        size_t existing_dim = existing_base.size();
        size_t matrix_dim = existing_dim + 1;
        if (matrix_dim > vector_space_dimension)
        {
            return false;
        }

        // Construct the basic matrix
        std::vector<std::vector<translation_coordinate_type>> base_matrix(existing_dim + 1, std::vector<translation_coordinate_type>(vector_space_dimension));
        for (size_t i = 0; i < existing_dim; i++)
        {
            for (size_t j = 0; j < vector_space_dimension; j++)
            {
                base_matrix[i][j] = existing_base[i].vec[j];
            }
        }

        std::vector<std::vector<translation_coordinate_type>> Mt_M(matrix_dim, std::vector<translation_coordinate_type>(matrix_dim, 0));

#pragma omp parallel for
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            for (size_t j = 0; j < vector_space_dimension; j++)
            {
                for (size_t k = 0; k < existing_dim; k++)
                {
                    Mt_M[i][j] += base_matrix[k][i] * base_matrix[k][j];
                }
            }
        }

        translation_coordinate_type MT_M_determinant = det(Mt_M);

        // If the absolute determinant result is below the threshold, then the vectors are linearly dependent
        if (translation_abs_function(MT_M_determinant) <= translation_coordinate_precision)
        {
            return false;
        }

        // Vectors are linearly independent
        return true;
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

    std::vector<ComponentInfo> PercolationGraph::get_component_percolation_info() const
    {
        // Obtain component decomposition
        std::vector<ComponentInfo> component_info = get_components();

        const int64_t comp_count = component_info.size();
        const size_t num_vertices = this->vertices.size();
        std::vector<bool> visited(num_vertices, false);
        std::vector<struct TranslationVector> original_positions(num_vertices);

        TranslationVector origin;
        for (size_t i = 0; i < vector_space_dimension; i++)
        {
            origin.vec[i] = 0;
        }

#pragma omp parallel for
        for (int64_t curr_component = 0; curr_component < comp_count; curr_component++)
        {
            // Determine percolation dimension for each of the

            ComponentInfo &curr_info = component_info[curr_component];
            size_t start_vertex = curr_info.vertices[0].index;

            std::vector<TranslationVector> basis_set;

            std::queue<std::pair<size_t, TranslationVector>, std::set<std::pair<size_t, TranslationVector>>> vertex_queue;
            vertex_queue.push({start_vertex, origin});

            while (!vertex_queue.empty())
            {
                auto vertex_position = vertex_queue.front();
                vertex_queue.pop();

                size_t vert_index = vertex_position.first;
                const TranslationVector &curr_position = vertex_position.second;

                // Possibly encountered different copy of vertex
                if (visited[vert_index])
                {
                    TranslationVector difference = curr_position - original_positions[vert_index];
                    // got new entry to basis set
                    if (check_translation_independent(basis_set, difference))
                    {
                        basis_set.push_back(difference);

                        // Maximum dimension attained, stop analysis
                        if (basis_set.size() >= vector_space_dimension)
                        {
                            break;
                        }
                    }
                    continue;
                }

                // Deal with first copy of vertex
                visited[vert_index] = true;
                original_positions[vert_index] = curr_position;
                for (auto edge : this->edges[vert_index])
                {
                    size_t neighbor = edge.first;
                    TranslationVector target_position = curr_position + edge.second.translation;
                    if (!visited[neighbor])
                    {
                        vertex_queue.push({neighbor, target_position});
                    }
                }
            }
            curr_info.perolation_dim = basis_set.size();
        }

        return component_info;
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