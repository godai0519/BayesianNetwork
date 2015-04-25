#ifndef BNI_EVALUATION_AIC_HPP
#define BNI_EVALUATION_AIC_HPP

#include <bayesian/graph.hpp>
#include <bayesian/evaluation/basic_info_criteria.hpp>

namespace bn {
namespace evaluation {
  
struct aic : basic_info_criteria {
    aic(sampler& sampling)
        : basic_info_criteria(sampling)
    {
    }

    double operator() (graph_t const& graph) const override
    {
        auto const likelihood = calc_likelihood(graph);
        auto const parameters = calc_parameters(graph);
        return likelihood + parameters;
    }
};

} // namespace evaluation
} // namespace bn

#endif // #ifndef BNI_EVALUATION_AIC_HPP
