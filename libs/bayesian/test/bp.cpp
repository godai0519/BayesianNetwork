#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"

// Pearl's Belief Propagation Algorithm -> Examples of application of algorithm
// http://www.cse.unsw.edu.au/~cs9417ml/Bayes/Pages/PearlPropagation.html
// 2 Samples
bn::graph_t make_pearl_bp_graph()
{
    bn::graph_t graph;
    auto vertex_r = graph.add_vertex();
    auto vertex_s = graph.add_vertex();
    auto vertex_w = graph.add_vertex();
    auto vertex_h = graph.add_vertex();
    auto edge_rw = graph.add_edge(vertex_r, vertex_w);
    auto edge_rh = graph.add_edge(vertex_r, vertex_h);
    auto edge_sh = graph.add_edge(vertex_s, vertex_h);
    
    {
        vertex_r->id = 1;
        vertex_r->selectable_num = 2;
        vertex_r->cpt.assign({}, vertex_r);

        bn::condition_t const cond;
        vertex_r->cpt[cond].second = {0.2, 0.8};
    }
    {
        vertex_s->id = 2;
        vertex_s->selectable_num = 2;
        vertex_s->cpt.assign({}, vertex_s);

        bn::condition_t const cond;
        vertex_s->cpt[cond].second = {0.1, 0.9};
    }
    {
        vertex_w->id = 3;
        vertex_w->selectable_num = 2;
        vertex_w->cpt.assign({vertex_r}, vertex_w);

        bn::condition_t const cond0 = {{vertex_r, 0}};
        bn::condition_t const cond1 = {{vertex_r, 1}};
        vertex_w->cpt[cond0].second = {1.0, 0.0};
        vertex_w->cpt[cond1].second = {0.2, 0.8};
    }
    {
        vertex_h->id = 4;
        vertex_h->selectable_num = 2;
        vertex_h->cpt.assign({vertex_r, vertex_s}, vertex_h);

        bn::condition_t const cond00 = {{vertex_r, 0}, {vertex_s, 0}};
        bn::condition_t const cond01 = {{vertex_r, 0}, {vertex_s, 1}};
        bn::condition_t const cond10 = {{vertex_r, 1}, {vertex_s, 0}};
        bn::condition_t const cond11 = {{vertex_r, 1}, {vertex_s, 1}};
        vertex_h->cpt[cond00].second = {1.0, 0.0};
        vertex_h->cpt[cond01].second = {0.9, 0.1};
        vertex_h->cpt[cond10].second = {1.0, 0.0};
        vertex_h->cpt[cond11].second = {0.0, 1.0};
    }

    return graph;
}

BOOST_AUTO_TEST_CASE( bp_pearl_part1 )
{
    bn::graph_t graph = make_pearl_bp_graph();
    auto const vertex = graph.vertex_list();
    std::vector<std::vector<double>> const teacher_r = {
        { 0.200, 0.800 },
        { 0.100, 0.900 },
        { 0.360, 0.640 },
        { 0.262, 0.738 }
    };

    bn::bp bp(graph);
    auto const result = bp();
    
    for(std::size_t i = 0; i < vertex.size(); ++i)
    {
        matrix_type const data = result.at(vertex[i]);

        BOOST_CHECK(data.height() == 1);

        for(std::size_t j = 0; j < data.width(); ++j)
        {
            BOOST_CHECK_CLOSE(data[0][j], teacher_r[i][j], 0.01);
        }
    }
}

BOOST_AUTO_TEST_CASE( bp_pearl_part2 )
{
    bn::graph_t graph = make_pearl_bp_graph();
    auto const vertex = graph.vertex_list();
    std::vector<std::vector<double>> const teacher_r = {
        { 0.7353, 0.2647 },
        { 0.3382, 0.6618 },
        { 0.7882, 0.2118 },
        { 1.0000, 0.0000 }
    };

    std::unordered_map<bn::vertex_type, bn::matrix_type> precondition;
    precondition[vertex[3]].resize(1, 2, 0);
    precondition[vertex[3]][0][0] = 1;

    bn::bp bp(graph);
    auto const result = bp(precondition);
    
    for(std::size_t i = 0; i < vertex.size(); ++i)
    {
        matrix_type const data = result.at(vertex[i]);

        BOOST_CHECK(data.height() == 1);

        for(std::size_t j = 0; j < data.width(); ++j)
        {
            BOOST_CHECK_CLOSE(data[0][j], teacher_r[i][j], 0.01);
        }
    }
}

//
// resume:
// 4 samples
//
/*
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

    edge_ab->likelihood = {true, mat_ba};
    edge_bc->likelihood = {true, mat_cb};
    edge_cd->likelihood = {true, mat_dc};
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = {true, mat_a};

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

    edge_ab->likelihood = {true, mat_ba};
    edge_bc->likelihood = {true, mat_cb};
    edge_cd->likelihood = {true, mat_dc};
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = {true, mat_a};

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

    edge_ab->likelihood = {true, mat_ba};
    edge_bc->likelihood = {true, mat_cb};
    edge_cd->likelihood = {true, mat_dc};
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = {true, mat_a};

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

    edge_ab->likelihood = {true, mat_ba};
    edge_bc->likelihood = {true, mat_cb};
    edge_cd->likelihood = {true, mat_dc};
    vertex_a->selectable_num = 3;
    vertex_b->selectable_num = 3;
    vertex_c->selectable_num = 2;
    vertex_d->selectable_num = 3;
    vertex_a->evidence = {true, mat_a};

    bn::bp func;
    auto const result = func(graph, vertex_a, {std::make_pair(vertex_d, 2)}); // P(A|D=3)
    for(int i = 0; i < 3; ++i)
    {
        // 誤差 3%検査
        BOOST_CHECK_CLOSE(result[i][0], mat_t[i][0], 3);
    }
}
*/
