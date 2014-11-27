#ifndef BNI_GRAPH_HPP
#define BNI_GRAPH_HPP

#include <vector>
#include <cstdint>
#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>
#include "matrix.hpp"

namespace bn {

struct vertex {
    int id;
    int selectable_num = 0; // æ‚è‚¤‚é’l‚Ì”
    boost::optional<matrix_type> evidence;
};

struct edge {
    boost::optional<matrix_type> likelihood;
};

struct graph_tag {
    std::string name;
};

// —vŒŸ“¢
typedef boost::adjacency_list<
    boost::listS,
    boost::vecS,
    boost::bidirectionalS,
    vertex,
    edge,
    graph_tag,
    boost::listS
> graph_t;

} // namespace bn

#endif // #ifndef BNI_MATRIX_HPP

