/**
* @file algorithm.hpp
* @brief Provide algorithms for networks.
* @author godai_0519
* @date 10/11/2017
*/

#ifndef BAYESIAN_NETWORKS_ALGORITHM_HPP
#define BAYESIAN_NETWORKS_ALGORITHM_HPP

#include <functional>
#include <vector>
#include <bayesian/network.hpp>
#include <bayesian/network/component.hpp>

namespace bn {

//! Generate a recursive function of the given function.
/*!
 *  @tparam      F: a type of the recursive function.
 *  @param[in]   f: a function to make the recursive function.
 *  @return      Generated recursive function.
**/
template<class F>
auto recursive(F f)
{
    return [f](auto... a) { return f(f, a...); };
}

//! Perform topological sort for given graph (however, a reversed series can be obtain).
/*! The topological sort gets a node series such that if a edge (u,v) is exist, v appears after u.
 *  This function is create a ``reversed'' topological sort series.
 *  The time complexity is O(V + E).
 *
 *  @tparam RepresentMethod: bn::adjacency_list or adjacency_matrix.
 *  @param[in]   g: a network including nodes which are targets of sort.
 *  @return      A node series created by ``reversed'' topological sort.
**/
template<class RepresentMethod>
std::vector<component::node_ptr> topological_sort(network<RepresentMethod> const& g)
{
    std::vector<component::node_ptr> result;
    result.reserve(g.all_node().size());

    // 0: no visit
    // 1: temporary
    // 2: permanent
    std::unordered_map<component::node_ptr, std::uint_fast8_t> marks;
    for(auto const& n : g.all_node())
        marks[n] = 0;

    // return true if a sort series of given node and its children is found;
    // return false if it is not found.
    auto visit = recursive(
        [&marks, &g, &result](auto visit, component::node_ptr const& n) -> bool
        {
            if(marks[n] == 1) return false; // the network is not DAG.
            else if(marks[n] == 0)
            {
                marks[n] = 1;
                for(auto const& child : g.child_nodes(n))
                    if(!visit(visit, child)) return false;
                marks[n] = 2;
                result.push_back(n);
            }
        });

    for(auto const& n : g.all_node())
        if(marks[n] == 0)
            if(!visit(n)) break;

    return result;
}

//! Judge that given graph is directed acyclic graph (DAG).
/*! Powered by toporogical sort, so the time complexity is O(V + E).
 *
 *  @tparam RepresentMethod: bn::adjacency_list or adjacency_matrix.
 *  @param[in]   g: a network to check that it is DAG.
 *  @return      true if the network is DAG.
**/
template<class RepresentMethod>
bool is_dag(network<RepresentMethod> const& g)
{
    return topological_sort(g).size() == g.all_node().size();
}

} // namespace bn

#endif // BAYESIAN_NETWORKS_ALGORITHM_HPP

