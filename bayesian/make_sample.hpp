/**
* @file make_sample.hpp
* @brief Make a random sample from network and cpts.
* @author godai_0519
* @date 10/11/2017
*/

#ifndef BAYESIAN_NETWORKS_MAKE_SAMPLE_HPP
#define BAYESIAN_NETWORKS_MAKE_SAMPLE_HPP

#include <random>
#include <vector>
#include <bayesian/network.hpp>
#include <bayesian/cpt.hpp>

namespace bn {

//! Perform topological sort for given graph (however, a reversed series can be obtain).
/*! The topological sort gets a node series such that if a edge (u,v) is exist, v appears after u.
 *  This function is create a ``reversed'' topological sort series.
 *  The time complexity is O(V + E).
 *
 *  @param[in]   g: a network including nodes which are targets of sort.
 *  @return      A node series created by ``reversed'' topological sort.
**/
template<class RepresentMethod>
cpt::condition_type make_sample(network<RepresentMethod> const& g, cpt_manager const& cpts, std::vector<component::node_ptr> const& sorted_nodes)
{
    // TODO: Policy
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());

    cpt::condition_type sample;
    for(auto it = sorted_nodes.crbegin(); it != sorted_nodes.crend(); ++it)
    {
        auto const probabilities = cpts.at(*it).at(sample);
        std::discrete_distribution<std::size_t> distribution(
            probabilities.cbegin(), probabilities.cend());
        sample[(*it)->get()] = distribution(engine);
    }

    return sample;
}

template<class RepresentMethod>
std::pair<cpt::condition_type, double> make_weighted_sample(network<RepresentMethod> const& g, cpt_manager const& cpts, std::unordered_map<bn::component::random_variable_ptr, std::size_t> const& evidence_nodes, std::vector<component::node_ptr> const& sorted_nodes)
{
    // TODO: Policy
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());

    cpt::condition_type sample;
    double weighted = 1.0;

    for(auto it = sorted_nodes.crbegin(); it != sorted_nodes.crend(); ++it)
    {
        if(evidence_nodes.find((*it)->get()) == evidence_nodes.end())
        {
            // it is NOT a evidence node.
            auto const probabilities = cpts.at(*it).at(sample);
            std::discrete_distribution<std::size_t> distribution(
                    probabilities.cbegin(), probabilities.cend());
            sample[(*it)->get()] = distribution(engine);
        }
        else
        {
            // it is a evidence node.
            weighted *= cpts.at(*it).at(sample)[evidence_nodes.at((*it)->get())];
            sample[(*it)->get()] = evidence_nodes.at((*it)->get());
        }
    }

    return { std::move(sample), weighted };
}

} // namespace bn

#endif // BAYESIAN_NETWORKS_MAKE_SAMPLE_HPP

