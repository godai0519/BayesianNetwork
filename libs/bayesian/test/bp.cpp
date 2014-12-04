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

    bn::matrix_type mat_a (1, 3);
    bn::matrix_type mat_ba(3, 3);
    bn::matrix_type mat_cb(3, 2);
    bn::matrix_type mat_dc(2, 3);
    bn::matrix_type mat_t (3, 1);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
    auto vertex_a = graph.add_vertex();
    auto vertex_b = graph.add_vertex();
    auto vertex_c = graph.add_vertex();
    auto vertex_d = graph.add_vertex();
    auto edge_ab = graph.add_edge(vertex_a, vertex_b);
    auto edge_bc = graph.add_edge(vertex_b, vertex_c);
    auto edge_cd = graph.add_edge(vertex_c, vertex_d);

    edge_ab->likelihood = mat_ba;
    edge_bc->likelihood = mat_cb;
    edge_cd->likelihood = mat_dc;
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = mat_a;

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

    bn::matrix_type mat_a (1, 3);
    bn::matrix_type mat_ba(3, 3);
    bn::matrix_type mat_cb(3, 2);
    bn::matrix_type mat_dc(2, 3);
    bn::matrix_type mat_t (3, 1);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
    auto vertex_a = graph.add_vertex();
    auto vertex_b = graph.add_vertex();
    auto vertex_c = graph.add_vertex();
    auto vertex_d = graph.add_vertex();
    auto edge_ab = graph.add_edge(vertex_a, vertex_b);
    auto edge_bc = graph.add_edge(vertex_b, vertex_c);
    auto edge_cd = graph.add_edge(vertex_c, vertex_d);

    edge_ab->likelihood = mat_ba;
    edge_bc->likelihood = mat_cb;
    edge_cd->likelihood = mat_dc;
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = mat_a;

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

    bn::matrix_type mat_a (1, 3);
    bn::matrix_type mat_ba(3, 3);
    bn::matrix_type mat_cb(3, 2);
    bn::matrix_type mat_dc(2, 3);
    bn::matrix_type mat_t (3, 1);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
    auto vertex_a = graph.add_vertex();
    auto vertex_b = graph.add_vertex();
    auto vertex_c = graph.add_vertex();
    auto vertex_d = graph.add_vertex();
    auto edge_ab = graph.add_edge(vertex_a, vertex_b);
    auto edge_bc = graph.add_edge(vertex_b, vertex_c);
    auto edge_cd = graph.add_edge(vertex_c, vertex_d);

    edge_ab->likelihood = mat_ba;
    edge_bc->likelihood = mat_cb;
    edge_cd->likelihood = mat_dc;
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = mat_a;

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

    bn::matrix_type mat_a (1, 3);
    bn::matrix_type mat_ba(3, 3);
    bn::matrix_type mat_cb(3, 2);
    bn::matrix_type mat_dc(2, 3);
    bn::matrix_type mat_t (3, 1);

    mat_a .assign(source_a .begin(), source_a .end());
    mat_ba.assign(source_ba.begin(), source_ba.end());
    mat_cb.assign(source_cb.begin(), source_cb.end());
    mat_dc.assign(source_dc.begin(), source_dc.end());
    mat_t .assign(  teacher.begin(),   teacher.end());

    bn::graph_t graph;
    auto vertex_a = graph.add_vertex();
    auto vertex_b = graph.add_vertex();
    auto vertex_c = graph.add_vertex();
    auto vertex_d = graph.add_vertex();
    auto edge_ab = graph.add_edge(vertex_a, vertex_b);
    auto edge_bc = graph.add_edge(vertex_b, vertex_c);
    auto edge_cd = graph.add_edge(vertex_c, vertex_d);

    edge_ab->likelihood = mat_ba;
    edge_bc->likelihood = mat_cb;
    edge_cd->likelihood = mat_dc;
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = mat_a;

    bn::bp func;
    auto const result = func(graph, vertex_a, {std::make_pair(vertex_d, 2)}); // P(A|D=3)
    for(int i = 0; i < 3; ++i)
    {
        // 誤差 3%検査
        BOOST_CHECK_CLOSE(result[i][0], mat_t[i][0], 3);
    }
}

