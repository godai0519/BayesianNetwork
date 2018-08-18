#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <bayesian/network.hpp>
#include <bayesian/network/adjacency_list.hpp>
#include <bayesian/network/adjacency_matrix.hpp>

BOOST_AUTO_TEST_CASE(list_artificial_network)
{
    bn::network<bn::adjacency_list> network;

    // make nodes
    auto const node_a = network.add_node();
    auto const node_b = network.add_node();
    auto const node_c = network.add_node();
    auto const node_d = network.add_node();

    // make arcs
    auto const edge_ab = network.add_arc(node_a, node_b);
    auto const edge_bc = network.add_arc(node_b, node_c);
    auto const edge_cd = network.add_arc(node_c, node_d);

    // validity of nodes and arcs
    BOOST_CHECK(node_a && node_b && node_c && node_d && edge_ab && edge_bc && edge_cd);
    BOOST_CHECK(network.all_node().size() == 4);
    BOOST_CHECK(network.all_arc().size() == 3);

    // test of source/target
    BOOST_CHECK(network.source(edge_ab) == node_a);
    BOOST_CHECK(network.target(edge_ab) == node_b);
    BOOST_CHECK(network.source(edge_bc) == node_b);
    BOOST_CHECK(network.target(edge_bc) == node_c);
    BOOST_CHECK(network.source(edge_cd) == node_c);
    BOOST_CHECK(network.target(edge_cd) == node_d);

    // test of is_connect
    BOOST_CHECK(network.is_connect(node_a, edge_ab));
    BOOST_CHECK(network.is_connect(node_b, edge_ab));
    BOOST_CHECK(network.is_connect(node_b, edge_bc));
    BOOST_CHECK(network.is_connect(node_c, edge_bc));
    BOOST_CHECK(network.is_connect(node_c, edge_cd));
    BOOST_CHECK(network.is_connect(node_d, edge_cd));

    // test of is_adjacent
    BOOST_CHECK(!network.is_adjacent(node_a, node_a));
    BOOST_CHECK( network.is_adjacent(node_a, node_b));
    BOOST_CHECK(!network.is_adjacent(node_a, node_c));
    BOOST_CHECK(!network.is_adjacent(node_a, node_d));
    BOOST_CHECK(!network.is_adjacent(node_b, node_a));
    BOOST_CHECK(!network.is_adjacent(node_b, node_b));
    BOOST_CHECK( network.is_adjacent(node_b, node_c));
    BOOST_CHECK(!network.is_adjacent(node_b, node_d));
    BOOST_CHECK(!network.is_adjacent(node_c, node_a));
    BOOST_CHECK(!network.is_adjacent(node_c, node_b));
    BOOST_CHECK(!network.is_adjacent(node_c, node_c));
    BOOST_CHECK( network.is_adjacent(node_c, node_d));
    BOOST_CHECK(!network.is_adjacent(node_d, node_a));
    BOOST_CHECK(!network.is_adjacent(node_d, node_b));
    BOOST_CHECK(!network.is_adjacent(node_d, node_c));
    BOOST_CHECK(!network.is_adjacent(node_d, node_d));

    // make a node and arcs
    auto const node_e = network.add_node();
    auto const edge_ac = network.add_arc(node_a, node_c);
    auto const edge_eb = network.add_arc(node_e, node_b);
    auto const edge_ec = network.add_arc(node_e, node_c);
    BOOST_CHECK(node_e && edge_ac && edge_eb && edge_ec);
    BOOST_CHECK(network.all_node().size() == 5);
    BOOST_CHECK(network.all_arc().size() == 6);

    // parent/child node
    auto const parent_c = network.parent_nodes(node_c);
    auto const child_c = network.child_nodes(node_c);
    BOOST_CHECK(parent_c.size() == 3);
    BOOST_CHECK(child_c.size() == 1);
    BOOST_CHECK(std::find(parent_c.cbegin(), parent_c.cend(), node_a) != parent_c.cend());
    BOOST_CHECK(std::find(parent_c.cbegin(), parent_c.cend(), node_b) != parent_c.cend());
    BOOST_CHECK(std::find(parent_c.cbegin(), parent_c.cend(), node_e) != parent_c.cend());
    BOOST_CHECK(std::find(child_c.cbegin(), child_c.cend(), node_d) != child_c.cend());

    // remove arcs
    network.remove_arc(node_e, node_b);
    network.remove_arc(edge_ec);
    BOOST_CHECK(network.all_node().size() == 5);
    BOOST_CHECK(network.all_arc().size() == 4);

    // remove an isolated node
    network.remove_node(node_e);
    BOOST_CHECK(network.all_node().size() == 4);
    BOOST_CHECK(network.all_arc().size() == 4);

    // remove a node
    network.remove_node(node_c);
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 1);

    // add an arc
    auto const edge_ad = network.add_arc(node_a, node_d);
    BOOST_CHECK(edge_ad);
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 2);

    // remove all arcs
    network.remove_all_arc();
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 0);

    // add an arc
    auto const edge_bd = network.add_arc(node_b, node_d);
    BOOST_CHECK(edge_bd);
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 1);

    // remove all nodes
    network.remove_all_node();
    BOOST_CHECK(network.all_node().size() == 0);
    BOOST_CHECK(network.all_arc().size() == 0);
}

BOOST_DATA_TEST_CASE(list_add_full_arc, boost::unit_test::data::xrange(10) ^ boost::unit_test::data::random(3, 50), i, node_num)
{
    bn::network<bn::adjacency_list> network;

    // make nodes
    std::vector<bn::component::node_ptr> nodes;
    for (std::size_t i = 0; i < node_num; ++i)
        nodes.push_back(network.add_node());

    // make arcs
    std::vector<bn::component::arc_ptr> arcs;
    for (std::size_t i = 0; i < node_num; ++i) for (std::size_t j = i + 1; j < node_num; ++j)
        arcs.push_back(network.add_arc(nodes.at(i), nodes.at(j)));

    // validity of nodes and arcs
    for (auto const& node : nodes) BOOST_CHECK(node);
    for (auto const& arc : arcs) BOOST_CHECK(arc);
    BOOST_CHECK(network.all_node().size() == node_num);
    BOOST_CHECK(network.all_arc().size() == node_num * (node_num - 1) / 2);

    // remove a node
    network.remove_node(nodes.back());
    BOOST_CHECK(network.all_node().size() == node_num - 1);
    BOOST_CHECK(network.all_arc().size() == node_num * (node_num - 1) / 2 - (node_num - 1));

    // remove a arc
    network.remove_arc(arcs.front());
    BOOST_CHECK(network.all_node().size() == node_num - 1);
    BOOST_CHECK(network.all_arc().size() == node_num * (node_num - 1) / 2 - (node_num - 1) - 1);

    // remove all nodes
    network.remove_all_node();
    BOOST_CHECK(network.all_node().size() == 0);
    BOOST_CHECK(network.all_arc().size() == 0);
}

BOOST_AUTO_TEST_CASE(list_change_direction)
{
    bn::network<bn::adjacency_list> network;
    std::vector<bn::component::arc_ptr> arcs;

    // Initialize
    auto const node_a = network.add_node();
    auto const node_b = network.add_node();
    auto const node_c = network.add_node();
    auto const node_d = network.add_node();
    auto const node_e = network.add_node();
    arcs.push_back(network.add_arc(node_a, node_b));
    arcs.push_back(network.add_arc(node_a, node_c));
    arcs.push_back(network.add_arc(node_b, node_c));
    arcs.push_back(network.add_arc(node_c, node_d));
    arcs.push_back(network.add_arc(node_e, node_b));
    arcs.push_back(network.add_arc(node_e, node_c));

    for(auto const& arc : arcs)
    {
        auto const old_source = network.source(arc);
        auto const old_target = network.target(arc);
        network.change_direction(arc);

        BOOST_CHECK(network.source(arc) == old_target);
        BOOST_CHECK(network.target(arc) == old_source);
    }
}

BOOST_AUTO_TEST_CASE(list_clone_node)
{
    bn::network<bn::adjacency_list> network;
    auto const node_a = network.add_node();
    auto const node_b = network.add_clone_node(node_a);
    node_a->get()->max_value = 5;

    BOOST_CHECK(node_a && node_b);
    BOOST_CHECK(node_b->get()->max_value == 5);
    BOOST_CHECK(node_a->get() == node_b->get());
    BOOST_CHECK(network.all_node().size() == 2);

    auto const edge_ab = network.add_arc(node_a, node_b);
    BOOST_CHECK(edge_ab);
    BOOST_CHECK(network.all_arc().size() == 1);
    BOOST_CHECK(network.source(edge_ab) == node_a);
    BOOST_CHECK(network.target(edge_ab) == node_b);

    network.remove_arc(edge_ab);
    BOOST_CHECK(network.all_node().size() == 2);
    BOOST_CHECK(network.all_arc().size() == 0);

    network.remove_node(node_a);
    BOOST_CHECK(network.all_node().size() == 1);
    BOOST_CHECK(network.all_arc().size() == 0);
}


BOOST_AUTO_TEST_CASE(list_move_clone)
{
    bn::network<bn::adjacency_list> network;
    auto const node_a = network.add_node();
    auto const node_b = network.add_node();
    auto const node_c = network.add_node();
    auto const edge_ab = network.add_arc(node_a, node_b);
    auto const edge_bc = network.add_arc(node_b, node_c);
    auto const edge_ca = network.add_arc(node_c, node_a);

    bn::network<bn::adjacency_list> cloned = network.clone();
    auto const cloned_nodes = cloned.all_node();
    auto const cloned_arcs = cloned.all_arc();
    BOOST_CHECK(std::find(cloned_nodes.cbegin(), cloned_nodes.cend(), node_a) != cloned_nodes.cend());
    BOOST_CHECK(std::find(cloned_nodes.cbegin(), cloned_nodes.cend(), node_b) != cloned_nodes.cend());
    BOOST_CHECK(std::find(cloned_nodes.cbegin(), cloned_nodes.cend(), node_c) != cloned_nodes.cend());
    BOOST_CHECK(std::find(cloned_arcs.cbegin(), cloned_arcs.cend(), edge_ab) != cloned_arcs.cend());
    BOOST_CHECK(std::find(cloned_arcs.cbegin(), cloned_arcs.cend(), edge_bc) != cloned_arcs.cend());
    BOOST_CHECK(std::find(cloned_arcs.cbegin(), cloned_arcs.cend(), edge_ca) != cloned_arcs.cend());

    bn::network<bn::adjacency_list> moved = std::move(network);
    auto const moved_nodes = moved.all_node();
    auto const moved_arcs = moved.all_arc();
    BOOST_CHECK(std::find(moved_nodes.cbegin(), moved_nodes.cend(), node_a) != moved_nodes.cend());
    BOOST_CHECK(std::find(moved_nodes.cbegin(), moved_nodes.cend(), node_b) != moved_nodes.cend());
    BOOST_CHECK(std::find(moved_nodes.cbegin(), moved_nodes.cend(), node_c) != moved_nodes.cend());
    BOOST_CHECK(std::find(moved_arcs.cbegin(), moved_arcs.cend(), edge_ab) != moved_arcs.cend());
    BOOST_CHECK(std::find(moved_arcs.cbegin(), moved_arcs.cend(), edge_bc) != moved_arcs.cend());
    BOOST_CHECK(std::find(moved_arcs.cbegin(), moved_arcs.cend(), edge_ca) != moved_arcs.cend());
}

BOOST_AUTO_TEST_CASE(matrix_artificial_network)
{
    bn::network<bn::adjacency_matrix> network;

    // make nodes
    auto const node_a = network.add_node();
    auto const node_b = network.add_node();
    auto const node_c = network.add_node();
    auto const node_d = network.add_node();

    // make arcs
    auto const edge_ab = network.add_arc(node_a, node_b);
    auto const edge_bc = network.add_arc(node_b, node_c);
    auto const edge_cd = network.add_arc(node_c, node_d);

    // validity of nodes and arcs
    BOOST_CHECK(node_a && node_b && node_c && node_d && edge_ab && edge_bc && edge_cd);
    BOOST_CHECK(network.all_node().size() == 4);
    BOOST_CHECK(network.all_arc().size() == 3);

    // test of source/target
    BOOST_CHECK(network.source(edge_ab) == node_a);
    BOOST_CHECK(network.target(edge_ab) == node_b);
    BOOST_CHECK(network.source(edge_bc) == node_b);
    BOOST_CHECK(network.target(edge_bc) == node_c);
    BOOST_CHECK(network.source(edge_cd) == node_c);
    BOOST_CHECK(network.target(edge_cd) == node_d);

    // test of is_connect
    BOOST_CHECK(network.is_connect(node_a, edge_ab));
    BOOST_CHECK(network.is_connect(node_b, edge_ab));
    BOOST_CHECK(network.is_connect(node_b, edge_bc));
    BOOST_CHECK(network.is_connect(node_c, edge_bc));
    BOOST_CHECK(network.is_connect(node_c, edge_cd));
    BOOST_CHECK(network.is_connect(node_d, edge_cd));

    // test of is_adjacent
    BOOST_CHECK(!network.is_adjacent(node_a, node_a));
    BOOST_CHECK(network.is_adjacent(node_a, node_b));
    BOOST_CHECK(!network.is_adjacent(node_a, node_c));
    BOOST_CHECK(!network.is_adjacent(node_a, node_d));
    BOOST_CHECK(!network.is_adjacent(node_b, node_a));
    BOOST_CHECK(!network.is_adjacent(node_b, node_b));
    BOOST_CHECK(network.is_adjacent(node_b, node_c));
    BOOST_CHECK(!network.is_adjacent(node_b, node_d));
    BOOST_CHECK(!network.is_adjacent(node_c, node_a));
    BOOST_CHECK(!network.is_adjacent(node_c, node_b));
    BOOST_CHECK(!network.is_adjacent(node_c, node_c));
    BOOST_CHECK(network.is_adjacent(node_c, node_d));
    BOOST_CHECK(!network.is_adjacent(node_d, node_a));
    BOOST_CHECK(!network.is_adjacent(node_d, node_b));
    BOOST_CHECK(!network.is_adjacent(node_d, node_c));
    BOOST_CHECK(!network.is_adjacent(node_d, node_d));

    // make a node and arcs
    auto const node_e = network.add_node();
    auto const edge_ac = network.add_arc(node_a, node_c);
    auto const edge_eb = network.add_arc(node_e, node_b);
    auto const edge_ec = network.add_arc(node_e, node_c);
    BOOST_CHECK(node_e && edge_ac && edge_eb && edge_ec);
    BOOST_CHECK(network.all_node().size() == 5);
    BOOST_CHECK(network.all_arc().size() == 6);

    // parent/child node
    auto const parent_c = network.parent_nodes(node_c);
    auto const child_c = network.child_nodes(node_c);
    BOOST_CHECK(parent_c.size() == 3);
    BOOST_CHECK(child_c.size() == 1);
    BOOST_CHECK(std::find(parent_c.cbegin(), parent_c.cend(), node_a) != parent_c.cend());
    BOOST_CHECK(std::find(parent_c.cbegin(), parent_c.cend(), node_b) != parent_c.cend());
    BOOST_CHECK(std::find(parent_c.cbegin(), parent_c.cend(), node_e) != parent_c.cend());
    BOOST_CHECK(std::find(child_c.cbegin(), child_c.cend(), node_d) != child_c.cend());

    // remove arcs
    network.remove_arc(node_e, node_b);
    network.remove_arc(edge_ec);
    BOOST_CHECK(network.all_node().size() == 5);
    BOOST_CHECK(network.all_arc().size() == 4);

    // remove an isolated node
    network.remove_node(node_e);
    BOOST_CHECK(network.all_node().size() == 4);
    BOOST_CHECK(network.all_arc().size() == 4);

    // remove a node
    network.remove_node(node_c);
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 1);

    // add an arc
    auto const edge_ad = network.add_arc(node_a, node_d);
    BOOST_CHECK(edge_ad);
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 2);

    // remove all arcs
    network.remove_all_arc();
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 0);

    // add an arc
    auto const edge_bd = network.add_arc(node_b, node_d);
    BOOST_CHECK(edge_bd);
    BOOST_CHECK(network.all_node().size() == 3);
    BOOST_CHECK(network.all_arc().size() == 1);

    // remove all nodes
    network.remove_all_node();
    BOOST_CHECK(network.all_node().size() == 0);
    BOOST_CHECK(network.all_arc().size() == 0);
}

BOOST_DATA_TEST_CASE(matrix_add_full_arc, boost::unit_test::data::xrange(10) ^ boost::unit_test::data::random(3, 50), i, node_num)
{
    bn::network<bn::adjacency_matrix> network;

    // make nodes
    std::vector<bn::component::node_ptr> nodes;
    for (std::size_t i = 0; i < node_num; ++i)
        nodes.push_back(network.add_node());

    // make arcs
    std::vector<bn::component::arc_ptr> arcs;
    for (std::size_t i = 0; i < node_num; ++i) for (std::size_t j = i + 1; j < node_num; ++j)
        arcs.push_back(network.add_arc(nodes.at(i), nodes.at(j)));

    // validity of nodes and arcs
    for (auto const& node : nodes) BOOST_CHECK(node);
    for (auto const& arc : arcs) BOOST_CHECK(arc);
    BOOST_CHECK(network.all_node().size() == node_num);
    BOOST_CHECK(network.all_arc().size() == node_num * (node_num - 1) / 2);

    // remove a node
    network.remove_node(nodes.back());
    BOOST_CHECK(network.all_node().size() == node_num - 1);
    BOOST_CHECK(network.all_arc().size() == node_num * (node_num - 1) / 2 - (node_num - 1));

    // remove a arc
    network.remove_arc(arcs.front());
    BOOST_CHECK(network.all_node().size() == node_num - 1);
    BOOST_CHECK(network.all_arc().size() == node_num * (node_num - 1) / 2 - (node_num - 1) - 1);

    // remove all nodes
    network.remove_all_node();
    BOOST_CHECK(network.all_node().size() == 0);
    BOOST_CHECK(network.all_arc().size() == 0);
}

BOOST_AUTO_TEST_CASE(matrix_change_direction)
{
    bn::network<bn::adjacency_matrix> network;
    std::vector<bn::component::arc_ptr> arcs;

    // Initialize
    auto const node_a = network.add_node();
    auto const node_b = network.add_node();
    auto const node_c = network.add_node();
    auto const node_d = network.add_node();
    auto const node_e = network.add_node();
    arcs.push_back(network.add_arc(node_a, node_b));
    arcs.push_back(network.add_arc(node_a, node_c));
    arcs.push_back(network.add_arc(node_b, node_c));
    arcs.push_back(network.add_arc(node_c, node_d));
    arcs.push_back(network.add_arc(node_e, node_b));
    arcs.push_back(network.add_arc(node_e, node_c));

    for (auto const& arc : arcs)
    {
        auto const old_source = network.source(arc);
        auto const old_target = network.target(arc);
        network.change_direction(arc);

        BOOST_CHECK(network.source(arc) == old_target);
        BOOST_CHECK(network.target(arc) == old_source);
    }
}

BOOST_AUTO_TEST_CASE(matrix_clone_node)
{
    bn::network<bn::adjacency_matrix> network;
    auto const node_a = network.add_node();
    auto const node_b = network.add_clone_node(node_a);
    node_a->get()->max_value = 5;

    BOOST_CHECK(node_a && node_b);
    BOOST_CHECK(node_b->get()->max_value == 5);
    BOOST_CHECK(node_a->get() == node_b->get());
    BOOST_CHECK(network.all_node().size() == 2);

    auto const edge_ab = network.add_arc(node_a, node_b);
    BOOST_CHECK(edge_ab);
    BOOST_CHECK(network.all_arc().size() == 1);
    BOOST_CHECK(network.source(edge_ab) == node_a);
    BOOST_CHECK(network.target(edge_ab) == node_b);

    network.remove_arc(edge_ab);
    BOOST_CHECK(network.all_node().size() == 2);
    BOOST_CHECK(network.all_arc().size() == 0);

    network.remove_node(node_a);
    BOOST_CHECK(network.all_node().size() == 1);
    BOOST_CHECK(network.all_arc().size() == 0);
}


BOOST_AUTO_TEST_CASE(matrix_move_clone)
{
    bn::network<bn::adjacency_matrix> network;
    auto const node_a = network.add_node();
    auto const node_b = network.add_node();
    auto const node_c = network.add_node();
    auto const edge_ab = network.add_arc(node_a, node_b);
    auto const edge_bc = network.add_arc(node_b, node_c);
    auto const edge_ca = network.add_arc(node_c, node_a);


    bn::network<bn::adjacency_matrix> cloned = network.clone();
    auto const cloned_nodes = cloned.all_node();
    auto const cloned_arcs = cloned.all_arc();
    BOOST_CHECK(std::find(cloned_nodes.cbegin(), cloned_nodes.cend(), node_a) != cloned_nodes.cend());
    BOOST_CHECK(std::find(cloned_nodes.cbegin(), cloned_nodes.cend(), node_b) != cloned_nodes.cend());
    BOOST_CHECK(std::find(cloned_nodes.cbegin(), cloned_nodes.cend(), node_c) != cloned_nodes.cend());
    BOOST_CHECK(std::find(cloned_arcs.cbegin(), cloned_arcs.cend(), edge_ab) != cloned_arcs.cend());
    BOOST_CHECK(std::find(cloned_arcs.cbegin(), cloned_arcs.cend(), edge_bc) != cloned_arcs.cend());
    BOOST_CHECK(std::find(cloned_arcs.cbegin(), cloned_arcs.cend(), edge_ca) != cloned_arcs.cend());

    bn::network<bn::adjacency_matrix> moved = std::move(network);
    auto const moved_nodes = moved.all_node();
    auto const moved_arcs = moved.all_arc();
    BOOST_CHECK(std::find(moved_nodes.cbegin(), moved_nodes.cend(), node_a) != moved_nodes.cend());
    BOOST_CHECK(std::find(moved_nodes.cbegin(), moved_nodes.cend(), node_b) != moved_nodes.cend());
    BOOST_CHECK(std::find(moved_nodes.cbegin(), moved_nodes.cend(), node_c) != moved_nodes.cend());
    BOOST_CHECK(std::find(moved_arcs.cbegin(), moved_arcs.cend(), edge_ab) != moved_arcs.cend());
    BOOST_CHECK(std::find(moved_arcs.cbegin(), moved_arcs.cend(), edge_bc) != moved_arcs.cend());
    BOOST_CHECK(std::find(moved_arcs.cbegin(), moved_arcs.cend(), edge_ca) != moved_arcs.cend());
}
