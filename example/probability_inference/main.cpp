#include <iostream>
#include <bayesian/bayesian_network.hpp>
#include <bayesian/bp.hpp>
#include <bayesian/likelihood_weighting.hpp>

bn::graph_t make_network()
{
    bn::graph_t graph;

    // Add Nodes. (4 nodes)
    graph.add_vertex();
    graph.add_vertex();
    graph.add_vertex();
    graph.add_vertex();

    // Get nodes list.
    // stored in the order that invoked "add_vertex()"
    // so, "returned value of add_vertex" == graph.vertex_list[i]
    std::vector<bn::vertex_type> const& vertex_list = graph.vertex_list();

    // Set selectable_num.
    // The number of values that can be took by the random variable(node).
    // ex. if "NodeA" can be {0, 1, 2, 3}, then you set it 4.
    vertex_list[0]->selectable_num = vertex_list[2]-selectable_num = 2;
    vertex_list[1]->selectable_num = vertex_list[3]-selectable_num = 6;

    // Add Edges.
    // if invoke "add_edge(A, B)", then make edge "A -> B".
    graph.add_edge(node_E, node_T);
    graph.add_edge(node_T, node_T_next);
    graph.add_edge(node_E, node_T_next);
    graph.add_edge(node_E_next, node_T_next);

    // !!! Make CPT by your Samples !!!
    //
    // if you use CSV data file(1 row = 1 sample, 1 column == 1 random variable), you can use bn::bayesian_network::load_data.
    //   ex1. bn::bayesian_network<void> bn;
    //        bn.load_data("samples.csv", graph.vertex_list());
    //        bn.load_cpt(graph);
    //   ex2. bn::bayesian_network<void> bn;
    //        bn.load_cpt_by_save_memory("samples.csv", graph.vertex_list(), graph);
    //
    // otherwise, repeat below code
    //   vertex_list[i]->cpt.assign({parent_nodes}, vertex_list[i]);
    //   bn::condition_t const cond = {{parent_node1, ?}, {parent_node2, ?}};
    //   vertex_list[i]->cpt[cond].second = {0.5, 0.3, ...};
    //

    // Set evidence node. In this example, vertex_list[i] is constantly 2. (and selectable_num is 4)
    std::unordered_map<bn::vertex_type, int> evidence_for_likelihood_weighting;
    evidence_for_likelihood_weighting[vertex_list[i]] = 2;

    std::unordered_map<bn::vertex_type, matrix_type> evidence_for_belief_propagation;
    evidence_for_belief_propagation[vertex_list[i]].resize(1, 4, 0.0);
    evidence_for_belief_propagation[vertex_list[i]][0][2] = 1;

    // Calculate! (Loopy Belief Propagation)
    bn::belief_propagation bp(graph);
    std::unordered_map<vertex_type, matrix_type> const result_bp =
        bp(evidence_for_belief_propagation, 0.000001/* Radius of convergence */);

    // Calculate! (Likelihood Weighting)
    bn::likelihood_weighting lw(graph);
    std::unordered_map<vertex_type, matrix_type> const result_bp =
        func_lw(evidence_for_likelihood_weighting, 1'000'000/* Num of auto sample */);
}

