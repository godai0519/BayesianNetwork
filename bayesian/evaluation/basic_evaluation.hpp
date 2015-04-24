#ifndef BNI_EVALUATION_BASIC_EVALUATION_HPP
#define BNI_EVALUATION_BASIC_EVALUATION_HPP

#include <bayesian/graph.hpp>

namespace bn {
namespace evaluation {
    
struct basic_evaluation {
    // API
    virtual double operator() (graph_t const& graph) const = 0;
};

} // namespace evaluation
} // namespace bn

#endif // #ifndef BNI_EVALUATION_BASIC_EVALUATION_HPP
