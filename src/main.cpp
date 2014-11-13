#include <iostream>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include "martix.hpp"

namespace bn {
    struct vertex {
        std::int_least32_t id;
        std::string        id_str; // twitterレスポンスに準じてみたがここに書くべきでない．辞書を作るべき．
    };

    struct edge {
        int cost;
        matrix likelihood;
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

    auto v1 = add_vertex(graph);
    auto v2 = add_vertex(graph);

    graph[v1].id = 1;
    graph[v2].id = 2;

    bool inserted = false;
    bn::graph_t::edge_descriptor e;
    boost::tie(e, inserted) = add_edge(v1, v2, graph);
    graph[e].cost = 100;

    std::vector<double> distance(boost::num_vertices(graph));
    boost::dijkstra_shortest_paths(
        graph, v1,
        boost::weight_map(boost::get(&bn::edge::cost, graph)).distance_map(&distance[0]));

    std::cout << "Tokyo-Nagoya : " << distance[v2] << "km" << std::endl;
}

