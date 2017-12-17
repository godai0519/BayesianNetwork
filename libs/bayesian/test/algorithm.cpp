#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <bayesian/network.hpp>
#include <bayesian/network/adjacency_list.hpp>
#include <bayesian/network/adjacency_matrix.hpp>
#include <bayesian/algorithm.hpp>

BOOST_AUTO_TEST_CASE(topological_sort_1)
{
    // v1 -> v2 -> v3 -> v4 -> v5
    bn::network<bn::adjacency_list> network;
    for(int i = 0; i < 5; ++i) network.add_node();
    for(int i = 0; i < 4; ++i)
        network.add_arc(
            network.all_node().at(i),
            network.all_node().at(i + 1));

    auto const topological = bn::topological_sort(network);
    auto const node = network.all_node();
    BOOST_CHECK(
        std::equal(
            topological.crbegin(), topological.crend(),
            node.cbegin(), node.cend()));

    BOOST_CHECK(bn::is_dag(network) == true);
}

BOOST_AUTO_TEST_CASE(topological_sort_2)
{
    bn::network<bn::adjacency_list> network;
    for(int i = 0; i < 6; ++i) network.add_node();
    network.add_arc(network.all_node().at(0), network.all_node().at(1));
    network.add_arc(network.all_node().at(0), network.all_node().at(2));
    network.add_arc(network.all_node().at(1), network.all_node().at(3));
    network.add_arc(network.all_node().at(2), network.all_node().at(3));
    network.add_arc(network.all_node().at(3), network.all_node().at(4));
    network.add_arc(network.all_node().at(3), network.all_node().at(5));

    auto const topological = bn::topological_sort(network);
    BOOST_CHECK(bn::is_dag(network) == true);
}

BOOST_AUTO_TEST_CASE(topological_sort_3)
{
    bn::network<bn::adjacency_list> network;
    for(int i = 0; i < 6; ++i) network.add_node();
    network.add_arc(network.all_node().at(0), network.all_node().at(2));
    network.add_arc(network.all_node().at(1), network.all_node().at(0));
    network.add_arc(network.all_node().at(2), network.all_node().at(3));
    network.add_arc(network.all_node().at(3), network.all_node().at(1));
    network.add_arc(network.all_node().at(3), network.all_node().at(4));
    network.add_arc(network.all_node().at(3), network.all_node().at(5));

    auto const topological = bn::topological_sort(network);
    BOOST_CHECK(bn::is_dag(network) == false);
}

