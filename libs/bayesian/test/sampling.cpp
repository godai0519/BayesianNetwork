#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/graph.hpp"
#include "bayesian/sampling.hpp"

BOOST_AUTO_TEST_CASE( sampling_standard )
{
    bn::graph_t graph;
    auto const vertex_1 = graph.add_vertex();
    auto const vertex_2 = graph.add_vertex();
    auto const vertex_3 = graph.add_vertex();
    auto const vertex_4 = graph.add_vertex();
    auto const vertex_5 = graph.add_vertex();
    auto edge_12 = graph.add_edge(vertex_1, vertex_2);
    auto edge_13 = graph.add_edge(vertex_1, vertex_3);
    auto edge_24 = graph.add_edge(vertex_2, vertex_4);
    auto edge_25 = graph.add_edge(vertex_2, vertex_5);
    auto edge_35 = graph.add_edge(vertex_3, vertex_5);

    {
        vertex_1->selectable_num = 2;
        vertex_1->cpt.assign({}, vertex_1);

        bn::condition_t const cond;
        vertex_1->cpt[cond].second = {0.5, 0.5};
    }
    {
        vertex_2->selectable_num = 2;
        vertex_2->cpt.assign({vertex_1}, vertex_2);

        bn::condition_t const cond_0 = {{vertex_1, 0}};
        bn::condition_t const cond_1 = {{vertex_1, 1}};
        vertex_2->cpt[cond_0].second = {0.8, 0.2};
        vertex_2->cpt[cond_1].second = {0.1, 0.9};
    }
    {
        vertex_3->selectable_num = 2;
        vertex_3->cpt.assign({vertex_1}, vertex_3);

        bn::condition_t const cond_0 = {{vertex_1, 0}};
        bn::condition_t const cond_1 = {{vertex_1, 1}};
        vertex_3->cpt[cond_0].second = {0.7, 0.3};
        vertex_3->cpt[cond_1].second = {0.4, 0.6};
    }
    {
        vertex_4->selectable_num = 2;
        vertex_4->cpt.assign({vertex_2}, vertex_4);

        bn::condition_t const cond_0 = {{vertex_2, 0}};
        bn::condition_t const cond_1 = {{vertex_2, 1}};
        vertex_4->cpt[cond_0].second = {0.6, 0.4};
        vertex_4->cpt[cond_1].second = {0.1, 0.9};
    }
    {
        vertex_5->selectable_num = 2;
        vertex_5->cpt.assign({vertex_2, vertex_3}, vertex_5);

        bn::condition_t const cond_00 = {{vertex_2, 0}, {vertex_3, 0}};
        bn::condition_t const cond_01 = {{vertex_2, 0}, {vertex_3, 1}};
        bn::condition_t const cond_10 = {{vertex_2, 1}, {vertex_3, 0}};
        bn::condition_t const cond_11 = {{vertex_2, 1}, {vertex_3, 1}};
        vertex_5->cpt[cond_00].second = {0.1, 0.9};
        vertex_5->cpt[cond_01].second = {0.2, 0.8};
        vertex_5->cpt[cond_10].second = {0.3, 0.7};
        vertex_5->cpt[cond_11].second = {0.4, 0.6};
    }

    bn::sampling func(graph);
    auto const result = func({{vertex_4,1}, {vertex_1, 0}});

    BOOST_CHECK_CLOSE(result.at(vertex_2)[0][0], 0.62, 10);
    BOOST_CHECK_CLOSE(result.at(vertex_2)[0][1], 0.38, 10);
}
