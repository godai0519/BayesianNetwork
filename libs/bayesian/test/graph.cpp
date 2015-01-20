#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/graph.hpp"

BOOST_AUTO_TEST_CASE( graph_trace_liner )
{
    bn::graph_t graph;
    auto vertex_a = graph.add_vertex();
    auto vertex_b = graph.add_vertex();
    auto vertex_c = graph.add_vertex();
    auto vertex_d = graph.add_vertex();
    auto edge_ab = graph.add_edge(vertex_a, vertex_b);
    auto edge_bc = graph.add_edge(vertex_b, vertex_c);
    auto edge_cd = graph.add_edge(vertex_c, vertex_d);

    BOOST_CHECK(graph.is_able_trace(vertex_a, vertex_d) == true);
    BOOST_CHECK(graph.is_able_trace(vertex_b, vertex_d) == true);
    BOOST_CHECK(graph.is_able_trace(vertex_c, vertex_d) == true);
    BOOST_CHECK(graph.is_able_trace(vertex_b, vertex_a) == false);
    BOOST_CHECK(graph.is_able_trace(vertex_c, vertex_a) == false);
    BOOST_CHECK(graph.is_able_trace(vertex_d, vertex_a) == false);
}

BOOST_AUTO_TEST_CASE( graph_trace_zigzag )
{
    bn::graph_t graph;
    auto vertex_a = graph.add_vertex();
    auto vertex_b = graph.add_vertex();
    auto vertex_c = graph.add_vertex();
    auto vertex_d = graph.add_vertex();
    auto edge_ab = graph.add_edge(vertex_a, vertex_b);
    auto edge_ac = graph.add_edge(vertex_a, vertex_c);
    auto edge_dc = graph.add_edge(vertex_d, vertex_c);

    BOOST_CHECK(graph.is_able_trace(vertex_a, vertex_d) == false);
    BOOST_CHECK(graph.is_able_trace(vertex_b, vertex_d) == false);
    BOOST_CHECK(graph.is_able_trace(vertex_d, vertex_b) == false);
}

BOOST_AUTO_TEST_CASE( graph_dag )
{
    bn::graph_t graph;
    auto vertex_a = graph.add_vertex();
    auto vertex_b = graph.add_vertex();
    auto vertex_c = graph.add_vertex();
    auto vertex_d = graph.add_vertex();

    // a -> b -> c -> d -> a
    BOOST_CHECK(graph.add_edge(vertex_a, vertex_b) != nullptr);
    BOOST_CHECK(graph.add_edge(vertex_b, vertex_c) != nullptr);
    BOOST_CHECK(graph.add_edge(vertex_c, vertex_d) != nullptr);
    BOOST_CHECK(graph.add_edge(vertex_d, vertex_a) == nullptr); // ERROR: nullptr
}
