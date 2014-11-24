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
        graph_t::vertex_descriptor const& target,
        std::vector<std::pair<graph_t::vertex_descriptor, int>> const& evidence
        );

private:
    matrix_type propagate_forward(
        graph_t const& graph,
        graph_t::vertex_descriptor const& target,
        std::vector<std::pair<graph_t::vertex_descriptor, int>> const& evidence
        );
    matrix_type propagate_backward(
        graph_t const& graph,
        graph_t::vertex_descriptor const& target,
        std::vector<std::pair<graph_t::vertex_descriptor, int>> const& evidence
        );
};

} // namespace bn

#endif // #ifndef BNI_BP_HPP

