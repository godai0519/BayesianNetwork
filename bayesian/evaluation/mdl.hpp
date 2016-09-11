#ifndef BNI_EVALUATION_MDL_HPP
#define BNI_EVALUATION_MDL_HPP

#include <bayesian/graph.hpp>
#include <bayesian/evaluation/basic_info_criteria.hpp>

namespace bn {
namespace evaluation {

struct mdl : basic_info_criteria {
    mdl(sampler const& sampling)
        : basic_info_criteria(sampling)
    {
    }

    double operator() (graph_t const& graph) const override
    {
        return (*this)(graph, graph.vertex_list());
    }

    double operator() (graph_t const& graph, std::vector<bn::vertex_type> const& vertex_list) const override
    {
        auto const likelihood = calc_likelihood(graph, vertex_list);
        auto const parameters = calc_parameters(graph);
        if(auto const sampling_size = this->sampling_size())
        {
            auto const correction = std::log2(sampling_size) / 2;
            return likelihood + parameters * correction;
        }
        else
        {
            throw std::runtime_error("Sampling is not finished yet.");
        }
    }
};

} // namespace evaluation
} // namespace bn

#endif // #ifndef BNI_EVALUATION_MDL_HPP
