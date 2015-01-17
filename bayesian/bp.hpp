#ifndef BNI_BP_HPP
#define BNI_BP_HPP

#include "graph.hpp"

namespace bn {

// 辺のlikelihoodを保持し，簡易に前の状態をコピーすることが可能にするため
class likelihood_list {
public:
    class value_type {
    public:
        value_type(vertex_type const& from, vertex_type const& to, edge_type const& edge);
        virtual ~value_type() = default;

        vertex_type const from() const;
        vertex_type const to() const;
        edge_type const edge() const;

        matrix_type& likelihood();
        matrix_type const& likelihood() const;

    private:
        vertex_type from_;
        vertex_type to_;
        edge_type edge_;
        matrix_type likelihood_;
    };

    matrix_type& add_manage(vertex_type const& from, vertex_type const& to, edge_type const& edge);
    void del_manage(edge_type const& edge);
    void del_manage(vertex_type const& from, vertex_type const& to);

    matrix_type& operator() (edge_type const& edge);
    matrix_type& operator() (vertex_type const& from, vertex_type const& to);
    matrix_type const& operator() (edge_type const& edge) const;
    matrix_type const& operator() (vertex_type const& from, vertex_type const& to) const;

private:
    std::vector<value_type>::iterator find(edge_type const& edge);
    std::vector<value_type>::iterator find(vertex_type const& from, vertex_type const& to);

    std::vector<value_type>::const_iterator find(edge_type const& edge) const;
    std::vector<value_type>::const_iterator find(vertex_type const& from, vertex_type const& to) const;

    std::vector<value_type> data_;
};

class bp {
public:
    bp() = default;
    virtual ~bp() = default;

    matrix_type operator()(
        graph_t const& graph,
        vertex_type const& target,
        std::vector<std::pair<vertex_type, int>> const& condition
        );

private:
    matrix_type propagate_forward(
        graph_t const& graph,
        vertex_type const& target,
        std::vector<std::pair<vertex_type, int>> const& condition
        );
    matrix_type propagate_backward(
        graph_t const& graph,
        vertex_type const& target,
        std::vector<std::pair<vertex_type, int>> const& condition
        );
};

} // namespace bn

#endif // #ifndef BNI_BP_HPP

