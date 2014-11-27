#ifndef BNI_BAYESIAN_NETWORK_HPP
#define BNI_BAYESIAN_NETWORK_HPP

#include "graph.hpp"

namespace bn {

template<class NodeType>
class bayesian_network {
public:
    template<class T>
    bayesian_network(T && graph);
    virtual ~bayesian_network() = default;

    template<class Func>
    auto apply(Func const& f) -> decltype(f(graph_));

private:
    graph_t graph_;
};

template<class NodeType> template<class T>
bayesian_network<NodeType>::bayesian_network(T && graph)
  : graph_(std::forward(graph))
{
}

template<class NodeType> template<class Func>
auto bayesian_network<NodeType>::apply(Func const& f) -> decltype(f(graph_))
{
    return f(graph);
}

} // namespace bn

#endif // #ifndef BNI_BAYESIAN_NETWORK_HPP

