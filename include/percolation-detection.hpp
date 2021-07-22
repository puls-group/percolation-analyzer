#include <vector>
#include <cmath>
#include <cstdlib>

// The value type of coordinates to be considered
using translation_coordinate_type = int64_t;
// If the absolute value of a translation coordinate entry is below this, it is considered to be zero
#define translation_coordinate_precision 0
// The function to obtain the absolute value for the coordinate type
#define translation_abs_function std::abs

// alternative values for floating point systems
//using translation_coordinate_type = float;
//#define translation_coordinate_precision 1e-7
//#define translation_abs_function std::abs

// The vector space dimension of the system
#define vector_space_dimension 3

namespace percolation
{

    struct TranslationVector
    {
        translation_coordinate_type vec[vector_space_dimension];

        struct TranslationVector operator-() const;
        struct TranslationVector operator+(struct TranslationVector const &other) const;
        struct TranslationVector operator-(struct TranslationVector const &other) const;
        bool operator==(struct TranslationVector const &other) const;
    };

    struct VertexData
    {
        size_t index;
    };

    struct EdgeData
    {
        EdgeData(const TranslationVector &trans) : translation(trans){}
        TranslationVector translation;

        struct EdgeData inverse() const;
    };

    struct ComponentInfo
    {
        size_t component_index;
        size_t perolation_dim;
        std::vector<VertexData> vertices;
    };
    
    /**
     * @brief Function to calculate the determinant of a square matrix
     * 
     * May not be super optimized but does the job. Only true implementation for dimensions up to 3
     * For higher dimensions, please provide own implementation in the source code or an error message will be thrown.
     * 
     * @param matrix The matrix whose determinant is supposed to be calculated
     * @return translation_coordinate_type The determinant value
     */
    translation_coordinate_type det(const std::vector<std::vector<translation_coordinate_type>> &matrix);

    /**
     * @brief Function to check if a given basis set of translation vectors and another translation vector are linearly independent
     * 
     * A zero vector will always result in a false result.
     * If a non-zero vector is passed as new_vector, the result will be true if the existing_base set is empty or if it is lineary independent of the existing_base.
     * Otherwise false will be the result.
     * 
     * If the result is true, adding the new_vector to existing_base will create a new basis set of vectors.
     * 
     * @param existing_base The pre-existing, linearly independent set of vectors or empty if none have been added yet
     * @param new_vector The vector to be checked for linear independence from the existing basis set
     * @return true The new_vector can be added to the set maintaining its basis properties
     * @return false The new_vector cannot be added to the set so that the result has vector basis properties.
     */
    bool check_translation_independent(const std::vector<TranslationVector> &existing_base, const TranslationVector &new_vector);

    class PercolationGraph
    {
    public:
        /**
         * @brief Reserve memory for the desired maximum number of vertices.
         * 
         * @param max_index The maximum index of a vertex (starting from zero) to be added.
         * @return true Memory allocation has been successful
         * @return false Memory allocation has failed
         */
        bool reserve_vertices(size_t max_index);

        /**
         * @brief Add vertex information to keep track of
         * 
         * @param vertex_index The index of the vertex to be annotated
         * @param vertex_data The data to be associated with the vertex
         * @return true The vertex data has been added successfully
         * @return false Something went wrong adding the vertex metadata
         */
        bool add_vertex(size_t vertex_index, const VertexData &vertex_data);

        /**
         * @brief Add an edge to the molecular graph including its translation vector
         * 
         * The translation vector is required to be pointing from the base to the head vertex, as direction of the translation is relevant to the overall grid structure.
         * Both the given translation as well as the inverse edge will be added, if possible.
         *  
         * @param vertex_index_base 
         * @param vertex_index_head 
         * @param edge_data 
         * @return true The edge has successfully been added.
         * @return false 
         */
        bool add_edge(size_t vertex_index_base, size_t vertex_index_head, const EdgeData &edge_data);

        /**
         * @brief Wrapper to directly provide the TranslationVector instead of an EdgeData object
         * 
         * Automatically wraps the TranslationVector in a EdgeData object before insertion
         * 
         * @param vertex_index_base 
         * @param vertex_index_head 
         * @param edge_trans 
         * @return true 
         * @return false 
         */
        bool add_edge(size_t vertex_index_base, size_t vertex_index_head, const TranslationVector &edge_trans);

        /**
         * @brief Get a list of all connected component of the current graph and their respective percolation information.
         * 
         * Will run a breadth first search to detect the components of the graph and then find the percolation dimension of 
         * each component by a more thorough bfs analyzing the grid structure of the graph
         * 
         * @return std::vector<ComponentInfo> 
         */
        std::vector<ComponentInfo> get_component_percolation_info() const;

    protected:
        /**
         * @brief Member to keep track of vertex information 
         */
        std::vector<VertexData> vertices;

        /**
         * @brief Member to keep track of the edges of the constructed graph
         */
        std::vector<std::vector<std::pair<size_t, EdgeData>>> edges;
        
        /**
         * @brief Get the connected components of the current graph.
         * 
         * @return std::vector<ComponentInfo> One component info entry for each detected component.
         */
        std::vector<ComponentInfo> get_components() const;
    };
}