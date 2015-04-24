#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/graph.hpp"
#include "bayesian/inference/belief_propagation.hpp"

// Pearl's Belief Propagation Algorithm -> Examples of application of algorithm
// http://www.cse.unsw.edu.au/~cs9417ml/Bayes/Pages/PearlPropagation.html
// 2 Samples
bn::graph_t make_pearl_belief_propagation_graph()
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
        vertex_h->cpt[cond01].second = {1.0, 0.0};
        vertex_h->cpt[cond10].second = {0.9, 0.1};
        vertex_h->cpt[cond11].second = {0.0, 1.0};
    }

    return graph;
}

BOOST_AUTO_TEST_CASE( belief_propagation_pearl_part1 )
{
    bn::graph_t graph = make_pearl_belief_propagation_graph();
    auto const vertex = graph.vertex_list();
    std::vector<std::vector<double>> const teacher_r = {
        { 0.200, 0.800 },
        { 0.100, 0.900 },
        { 0.360, 0.640 },
        { 0.272, 0.728 }
    };

    bn::inference::belief_propagation bp(graph);
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

// 多分なんか違う
BOOST_AUTO_TEST_CASE( belief_propagation_pearl_part2 )
{
    bn::graph_t graph = make_pearl_belief_propagation_graph();
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

    bn::inference::belief_propagation bp(graph);
    auto const result = bp(precondition);

    for(std::size_t i = 0; i < vertex.size(); ++i)
    {
        matrix_type const data = result.at(vertex[i]);

        BOOST_CHECK(data.height() == 1);

        for(std::size_t j = 0; j < data.width(); ++j)
        {
            BOOST_CHECK_CLOSE(data[0][j], teacher_r[i][j], 0.1);
        }
    }
}

//
// resume14:
// 5 samples
//
bn::graph_t make_resume_graph()
{
    bn::graph_t graph;
    auto const vertex_A = graph.add_vertex();
    auto const vertex_B = graph.add_vertex();
    auto const vertex_C = graph.add_vertex();
    auto const vertex_D = graph.add_vertex();
    auto const edge_AB = graph.add_edge(vertex_A, vertex_B);
    auto const edge_BC = graph.add_edge(vertex_B, vertex_C);
    auto const edge_CD = graph.add_edge(vertex_C, vertex_D);

    {
        bn::condition_t const cond;

        vertex_A->id = 1;
        vertex_A->selectable_num = 3;
        vertex_A->cpt.assign({}, vertex_A);
        vertex_A->cpt[cond].second = { 0.30, 0.60, 0.10 };
    }
    {
        bn::condition_t const cond0 = {{vertex_A, 0}};
        bn::condition_t const cond1 = {{vertex_A, 1}};
        bn::condition_t const cond2 = {{vertex_A, 2}};

        vertex_B->id = 2;
        vertex_B->selectable_num = 3;
        vertex_B->cpt.assign({vertex_A}, vertex_B);
        vertex_B->cpt[cond0].second = { 0.20, 0.30, 0.50 };
        vertex_B->cpt[cond1].second = { 0.30, 0.30, 0.40 };
        vertex_B->cpt[cond2].second = { 0.80, 0.10, 0.10 };
    }
    {
        bn::condition_t const cond0 = {{vertex_B, 0}};
        bn::condition_t const cond1 = {{vertex_B, 1}};
        bn::condition_t const cond2 = {{vertex_B, 2}};

        vertex_C->id = 3;
        vertex_C->selectable_num = 2;
        vertex_C->cpt.assign({vertex_B}, vertex_C);
        vertex_C->cpt[cond0].second = { 0.50, 0.50 };
        vertex_C->cpt[cond1].second = { 0.70, 0.30 };
        vertex_C->cpt[cond2].second = { 0.40, 0.60 };
    }
    {
        bn::condition_t const cond0 = {{vertex_C, 0}};
        bn::condition_t const cond1 = {{vertex_C, 1}};

        vertex_D->id = 4;
        vertex_D->selectable_num = 3;
        vertex_D->cpt.assign({vertex_C}, vertex_D);
        vertex_D->cpt[cond0].second = { 0.40, 0.30, 0.30 };
        vertex_D->cpt[cond1].second = { 0.20, 0.60, 0.20 };
    }

    return graph;
}

BOOST_AUTO_TEST_CASE( belief_propagation_resume_ex )
{
    std::vector<double> const teacher_C = { 0.570, 0.430 };

    bn::graph_t graph = make_resume_graph();
    auto const vertex = graph.vertex_list();

    std::unordered_map<bn::vertex_type, bn::matrix_type> precondition;
    precondition[vertex[1]].resize(1, 3);
    precondition[vertex[3]].resize(1, 3);
    precondition[vertex[1]][0] = { 0.0, 0.0, 1.0 };
    precondition[vertex[3]][0] = { 1.0, 0.0, 0.0 };

    bn::inference::belief_propagation func(graph);
    auto const result = func(precondition);

    auto const data = result.at(vertex[2]);
    BOOST_CHECK(data.height() == 1);

    for(std::size_t i = 0; i < data.width(); ++i)
    {
        BOOST_CHECK_CLOSE(data[0][i], teacher_C[i], 3);
    }
}

BOOST_AUTO_TEST_CASE( belief_propagation_resume_sample1 )
{
    std::vector<double> const teacher_B = { 0.330, 0.170, 0.500 };

    bn::graph_t graph = make_resume_graph();
    auto const vertex = graph.vertex_list();

    std::unordered_map<bn::vertex_type, bn::matrix_type> precondition;
    precondition[vertex[2]].resize(1, 2);
    precondition[vertex[2]][0] = { 0.0, 1.0 };

    bn::inference::belief_propagation func(graph);
    auto const result = func(precondition);

    auto const data = result.at(vertex[1]);
    BOOST_CHECK(data.height() == 1);

    for(std::size_t i = 0; i < data.width(); ++i)
    {
        BOOST_CHECK_CLOSE(data[0][i], teacher_B[i], 3);
    }
}

BOOST_AUTO_TEST_CASE( belief_propagation_resume_sample2 )
{
    std::vector<double> const teacher_B = { 0.310, 0.190, 0.500 };

    bn::graph_t graph = make_resume_graph();
    auto const vertex = graph.vertex_list();

    std::unordered_map<bn::vertex_type, bn::matrix_type> precondition;
    precondition[vertex[0]].resize(1, 3);
    precondition[vertex[2]].resize(1, 2);
    precondition[vertex[0]][0] = { 0.0, 1.0, 0.0 };
    precondition[vertex[2]][0] = { 0.0, 1.0 };

    bn::inference::belief_propagation func(graph);
    auto const result = func(precondition);

    auto const data = result.at(vertex[1]);
    BOOST_CHECK(data.height() == 1);

    for(std::size_t i = 0; i < data.width(); ++i)
    {
        BOOST_CHECK_CLOSE(data[0][i], teacher_B[i], 3);
    }
}

BOOST_AUTO_TEST_CASE( belief_propagation_resume_sample3 )
{
    std::vector<double> const teacher_A = { 0.300, 0.600, 0.100 };

    bn::graph_t graph = make_resume_graph();
    auto const vertex = graph.vertex_list();

    std::unordered_map<bn::vertex_type, bn::matrix_type> precondition;
    precondition[vertex[3]].resize(1, 3);
    precondition[vertex[3]][0] = { 0.0, 0.0, 1.0 };

    bn::inference::belief_propagation func(graph);
    auto const result = func(precondition);

    auto const data = result.at(vertex[0]);
    BOOST_CHECK(data.height() == 1);

    for(std::size_t i = 0; i < data.width(); ++i)
    {
        BOOST_CHECK_CLOSE(data[0][i], teacher_A[i], 3);
    }
}

BOOST_AUTO_TEST_CASE( belief_propagation_resume_sample4 )
{
    std::vector<double> const teacher_B = { 0.200, 0.300, 0.500 };

    bn::graph_t graph = make_resume_graph();
    auto const vertex = graph.vertex_list();

    std::unordered_map<bn::vertex_type, bn::matrix_type> precondition;
    precondition[vertex[0]].resize(1, 3);
    precondition[vertex[0]][0] = { 1.0, 0.0, 0.0 };

    bn::inference::belief_propagation func(graph);
    auto const result = func(precondition);

    auto const data = result.at(vertex[1]);
    BOOST_CHECK(data.height() == 1);

    for(std::size_t i = 0; i < data.width(); ++i)
    {
        BOOST_CHECK_CLOSE(data[0][i], teacher_B[i], 3);
    }
}
