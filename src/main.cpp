#include <iostream>
#include <string>
#include <cstdint>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/optional.hpp>
#include "matrix.hpp"

namespace bn {

struct vertex {
    std::int_least32_t id;
    std::string        id_str; // twitterレスポンスに準じてみたがここに書くべきでない．辞書を作るべき．

    bool evidence = false;
    std::int_least32_t evidence_select = 0;
};

struct edge {
    boost::optional<matrix_type> likelihood;
};

struct graph_tag {
    std::string name;
};

// 要検討
typedef boost::adjacency_list<
    boost::listS,
    boost::vecS,
    boost::directedS,
    vertex,
    edge,
    graph_tag,
    boost::listS
> graph_t;

} // namespace bn

int main()
{
    bn::graph_t graph;
    graph[boost::graph_bundle].name = "Bayesian Network";
}

