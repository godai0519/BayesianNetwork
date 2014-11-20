#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"

BOOST_AUTO_TEST_CASE( bp_specify_both_ends )
{
    std::vector<double> const source_ba = { 0.20, 0.30, 0.50, 0.30, 0.30, 0.40, 0.80, 0.10, 0.10};
    std::vector<double> const source_cb = { 0.50, 0.50, 0.70, 0.30, 0.40, 0.60};
    std::vector<double> const source_dc = { 0.40, 0.30, 0.30, 0.20, 0.60, 0.20};
    std::vector<double> const teacher   = { 0.31, 0.19, 0.50};

    bn::matrix_type mat_ba(boost::extents[3][3]);
    bn::matrix_type mat_cb(boost::extents[3][2]);
    bn::matrix_type mat_dc(boost::extents[2][3]);
    bn::matrix_type mat_t (boost::extents[3][1]);

    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    ke::graph_t graph;
    auto vertex_a = add_vertex(graph);
    auto vertex_b = add_vertex(graph);
    auto vertex_c = add_vertex(graph);
    auto vertex_d = add_vertex(graph);
    auto edge_ab = add_edge(vertex_a, vertex_b, graph);
    auto edge_bc = add_edge(vertex_b, vertex_c, graph);
    auto edge_cd = add_edge(vertex_c, vertex_d, graph);

    graph[edge_ab.first].likelihood = mat_ba;
    graph[edge_bc.first].likelihood = mat_cb;
    graph[edge_cd.first].likelihood = mat_dc;

    ke::bp func;
    auto const result = func(graph, {std::make_pair(vertex_a, 2), std::make_pair(vertex_c, 2)});
    for(int i = 0; i < 3; ++i)
    {
        // 誤差 0.00001%検査
        BOOST_CHECK_CLOSE(result[i][0], mat_t[i][0], 0.00001);
    }
}

