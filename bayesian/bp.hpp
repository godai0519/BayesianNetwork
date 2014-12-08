#ifndef BNI_BP_HPP
#define BNI_BP_HPP

#include "graph.hpp"

namespace bn {

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

