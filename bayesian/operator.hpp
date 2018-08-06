/**
* @file operator.hpp
* @brief Provide operator functions for networks.
* @author godai_0519
* @date 08/02/2018
*/

#ifndef BAYESIAN_NETWORKS_OPERATOR_HPP
#define BAYESIAN_NETWORKS_OPERATOR_HPP

#include <memory>
#include <vector>
#include <bayesian/network/component.hpp>

namespace bn {

template<class Network>
auto all_node(Network const& network)
{
    return network.all_node();
}

template<class Network>
bool is_contained(component::random_variable_ptr const& rv, Network const& network)
{
    auto const& nodes = all_node(network);
    return std::any_of(std::cbegin(nodes), std::cend(nodes),
        [&rv](auto const& n)
        {
            return n->get() == rv;
        });
}

template<class Network>
bool is_contained(component::node_ptr const& node, Network const& network)
{
    auto const& nodes = all_node(network);
    return std::any_of(std::cbegin(nodes), std::cend(nodes),
        [&node](auto const& n)
        {
            return n == node;
        });
}

template<class Elems, class Network>
bool is_contained(Elems const& elements, Network const& network)
{
    return std::all_of(std::cbegin(elements), std::cend(elements),
        [&network](auto const& element)
        {
            return is_contained(element, network);
        });
}

template<class Condition, class Evidence>
bool is_consistent(Condition const& condition, Evidence const& evidence)
{
    for(auto const& e : evidence)
    {
        if(condition.at(e.first) != e.second)
            return false; // inconsistent
    }

    return true;
}

} // namespace bn

#endif // BAYESIAN_NETWORKS_OPERATOR_HPP
