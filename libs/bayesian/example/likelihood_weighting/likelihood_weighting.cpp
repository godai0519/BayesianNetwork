/**
* @file rejection_sampling.cpp
* @brief An example which infers using Rejection Sampling.
* @author godai_0519
* @date 02/26/2018
*/

//
// In this example, we reason probabilities shown in the following webpage:
// A Brief Introduction to Graphical Models and Bayesian Networks (By Kevin Murphy, 1998)
// https://www.cs.ubc.ca/~murphyk/Bayes/bnintro.html
//
// Thank you, Dr. Murphy.
//

#include <iostream>
#include <bayesian/network.hpp>
#include <bayesian/network/adjacency_matrix.hpp>
#include <bayesian/cpt.hpp>
#include <bayesian/inference/likelihood_weighting.hpp>

int main()
{
    // =============================
    //     Making a network
    // =============================

    bn::network<bn::adjacency_matrix> network;
    auto const sprinkler = network.add_node();
    auto const wet_grass = network.add_node();
    auto const cloudy = network.add_node();
    auto const rain = network.add_node();

    network.add_arc(cloudy, sprinkler);
    network.add_arc(cloudy, rain);
    network.add_arc(sprinkler, wet_grass);
    network.add_arc(rain, wet_grass);

    sprinkler->get()->max_value
        = wet_grass->get()->max_value
        = cloudy->get()->max_value
        = rain->get()->max_value
        = 2;

    // =============================
    //     Making a CPTs
    // =============================

    bn::cpt_manager cpts;
    bn::cpt cpt_sprinkler(sprinkler->get(), {cloudy->get()}); // TODO: able to remove ``.get()''
    bn::cpt cpt_wet_grass(wet_grass->get(), {sprinkler->get(), rain->get()});
    bn::cpt cpt_cloudy(cloudy->get());
    bn::cpt cpt_rain(rain->get(), {cloudy->get()});

    cpt_sprinkler[{{cloudy->get(), 0}}] = std::vector<double>{0.5, 0.5};
    cpt_sprinkler[{{cloudy->get(), 1}}] = std::vector<double>{0.9, 0.1};
    cpt_wet_grass[{{sprinkler->get(), 0}, {rain->get(), 0}}] = std::vector<double>{1.0, 0.0};
    cpt_wet_grass[{{sprinkler->get(), 1}, {rain->get(), 0}}] = std::vector<double>{0.1, 0.9};
    cpt_wet_grass[{{sprinkler->get(), 0}, {rain->get(), 1}}] = std::vector<double>{0.1, 0.9};
    cpt_wet_grass[{{sprinkler->get(), 1}, {rain->get(), 1}}] = std::vector<double>{0.01, 0.99};
    cpt_cloudy[{}] = std::vector<double>{0.5, 0.5};
    cpt_rain[{{cloudy->get(), 0}}] = std::vector<double>{0.8, 0.2};
    cpt_rain[{{cloudy->get(), 1}}] = std::vector<double>{0.2, 0.8};

    cpts.enroll(sprinkler, cpt_sprinkler);
    cpts.enroll(wet_grass, cpt_wet_grass);
    cpts.enroll(cloudy, cpt_cloudy);
    cpts.enroll(rain, cpt_rain);

    // ==============================================
    //     Infering Pr(Sprinkler = 1 | WetGrass = 1)
    // ==============================================

    bn::inference::likelihood_weighting<bn::adjacency_matrix> rs1(network, cpts);
    rs1.set_query(std::vector<bn::component::random_variable_ptr>{sprinkler->get()})
       .set_evidence(std::unordered_map<bn::component::random_variable_ptr, std::size_t>{ {wet_grass->get(), 1}});

    rs1.run(1'000'000);
    auto const target_probability1 =
        rs1.probability().at(sprinkler->get())[std::vector<std::size_t>{1}]; // Pr(Sprinkler = 1 | WetGrass = 1)
    std::cout << target_probability1 << std::endl; // 0.429....

    rs1.run(1'000'000);
    auto const target_probability1_more =
        rs1.probability().at(sprinkler->get())[std::vector<std::size_t>{1}]; // more accurate Pr(Sprinkler = 1 | WetGrass = 1)
    std::cout << target_probability1_more << std::endl; // 0.429....


    // ==============================================
    //     Infering Pr(Rain = 1 | WetGrass = 1)
    // ==============================================

    bn::inference::likelihood_weighting<bn::adjacency_matrix> rs2(network, cpts);
    rs2.set_query(std::vector<bn::component::random_variable_ptr>{rain->get()})
       .set_evidence(std::unordered_map<bn::component::random_variable_ptr, std::size_t>{ {wet_grass->get(), 1}});

    rs2.run(1'000'000);
    auto const target_probability2 =
        rs2.probability().at(rain->get())[std::vector<std::size_t>{1}]; // Pr(Rain = 1 | WetGrass = 1)
    std::cout << target_probability2 << std::endl; // 0.708....
}

