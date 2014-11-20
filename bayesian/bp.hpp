#ifndef BNI_BP_HPP
#define BNI_BP_HPP

#include "graph.hpp"

namespace bn {

class bp {
public:
    bp() = default;
    virtual ~bp() = default;

    double operator()(graph_t const graph);
};

} // namespace bn

#endif // #ifndef BNI_BP_HPP

