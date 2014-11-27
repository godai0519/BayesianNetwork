#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"

BOOST_AUTO_TEST_CASE( bp_specify_both_ends )
{
    std::vector<double> const source_a  = { 0.30, 0.60, 0.10};
    std::vector<double> const source_ba = { 0.20, 0.30, 0.50, 0.30, 0.30, 0.40, 0.80, 0.10, 0.10};
    std::vector<double> const source_cb = { 0.50, 0.50, 0.70, 0.30, 0.40, 0.60};
    std::vector<double> const source_dc = { 0.40, 0.30, 0.30, 0.20, 0.60, 0.20};
    std::vector<double> const teacher   = { 0.3125, 0.1875, 0.5000};

    bn::matrix_type mat_a (boost::extents[1][3]);
    bn::matrix_type mat_ba(boost::extents[3][3]);
    bn::matrix_type mat_cb(boost::extents[3][2]);
    bn::matrix_type mat_dc(boost::extents[2][3]);
    bn::matrix_type mat_t (boost::extents[3][1]);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
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
    graph[vertex_a].selectable_num = 3;
    graph[vertex_b].selectable_num = 3;
    graph[vertex_c].selectable_num = 2;
    graph[vertex_d].selectable_num = 3;
    graph[vertex_a].evidence = mat_a;

    bn::bp func;
    auto const result = func(graph, vertex_b, {std::make_pair(vertex_a, 1), std::make_pair(vertex_c, 1)}); // P(B|A=2,C=2)
    for(int i = 0; i < 3; ++i)
    {
        // 誤差 2%検査
        BOOST_CHECK_CLOSE(result[i][0], mat_t[i][0], 2);
    }
}

BOOST_AUTO_TEST_CASE( bp_specify_only_upstream )
{
    std::vector<double> const source_a  = { 0.30, 0.60, 0.10};
    std::vector<double> const source_ba = { 0.20, 0.30, 0.50, 0.30, 0.30, 0.40, 0.80, 0.10, 0.10};
    std::vector<double> const source_cb = { 0.50, 0.50, 0.70, 0.30, 0.40, 0.60};
    std::vector<double> const source_dc = { 0.40, 0.30, 0.30, 0.20, 0.60, 0.20};
    std::vector<double> const teacher   = { 0.20, 0.30, 0.50};

    bn::matrix_type mat_a (boost::extents[1][3]);
    bn::matrix_type mat_ba(boost::extents[3][3]);
    bn::matrix_type mat_cb(boost::extents[3][2]);
    bn::matrix_type mat_dc(boost::extents[2][3]);
    bn::matrix_type mat_t (boost::extents[3][1]);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
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
    graph[vertex_a].selectable_num = 3;
    graph[vertex_b].selectable_num = 3;
    graph[vertex_c].selectable_num = 2;
    graph[vertex_d].selectable_num = 3;
    graph[vertex_a].evidence = mat_a;

    bn::bp func;
    auto const result = func(graph, vertex_b, {std::make_pair(vertex_a, 0)}); // P(B|A=1)
    for(int i = 0; i < 3; ++i)
    {
        // 誤差 2%検査
        BOOST_CHECK_CLOSE(result[i][0], mat_t[i][0], 2);
    }
}

BOOST_AUTO_TEST_CASE( bp_specify_only_downstream )
{
    std::vector<double> const source_a  = { 0.30, 0.60, 0.10};
    std::vector<double> const source_ba = { 0.20, 0.30, 0.50, 0.30, 0.30, 0.40, 0.80, 0.10, 0.10};
    std::vector<double> const source_cb = { 0.50, 0.50, 0.70, 0.30, 0.40, 0.60};
    std::vector<double> const source_dc = { 0.40, 0.30, 0.30, 0.20, 0.60, 0.20};
    std::vector<double> const teacher   = { 0.33, 0.17, 0.50};

    bn::matrix_type mat_a (boost::extents[1][3]);
    bn::matrix_type mat_ba(boost::extents[3][3]);
    bn::matrix_type mat_cb(boost::extents[3][2]);
    bn::matrix_type mat_dc(boost::extents[2][3]);
    bn::matrix_type mat_t (boost::extents[3][1]);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
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
    graph[vertex_a].selectable_num = 3;
    graph[vertex_b].selectable_num = 3;
    graph[vertex_c].selectable_num = 2;
    graph[vertex_d].selectable_num = 3;
    graph[vertex_a].evidence = mat_a;

    bn::bp func;
    auto const result = func(graph, vertex_b, {std::make_pair(vertex_c, 1)}); // P(B|C=2)
    for(int i = 0; i < 3; ++i)
    {
        // 誤差 3%検査
        BOOST_CHECK_CLOSE(result[i][0], mat_t[i][0], 3);
    }
}

BOOST_AUTO_TEST_CASE( bp_highest_node )
{
    std::vector<double> const source_a  = { 0.30, 0.60, 0.10};
    std::vector<double> const source_ba = { 0.20, 0.30, 0.50, 0.30, 0.30, 0.40, 0.80, 0.10, 0.10};
    std::vector<double> const source_cb = { 0.50, 0.50, 0.70, 0.30, 0.40, 0.60};
    std::vector<double> const source_dc = { 0.40, 0.30, 0.30, 0.20, 0.60, 0.20};
    std::vector<double> const teacher   = { 0.30, 0.60, 0.10};

    bn::matrix_type mat_a (boost::extents[1][3]);
    bn::matrix_type mat_ba(boost::extents[3][3]);
    bn::matrix_type mat_cb(boost::extents[3][2]);
    bn::matrix_type mat_dc(boost::extents[2][3]);
    bn::matrix_type mat_t (boost::extents[3][1]);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
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
    graph[vertex_a].selectable_num = 3;
    graph[vertex_b].selectable_num = 3;
    graph[vertex_c].selectable_num = 2;
    graph[vertex_d].selectable_num = 3;
    graph[vertex_a].evidence = mat_a;

    bn::bp func;
    auto const result = func(graph, vertex_a, {std::make_pair(vertex_d, 2)}); // P(A|D=3)
    for(int i = 0; i < 3; ++i)
    {
        // 誤差 3%検査
        BOOST_CHECK_CLOSE(result[i][0], mat_t[i][0], 3);
    }
}

