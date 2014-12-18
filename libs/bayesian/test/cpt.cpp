#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/graph.hpp"

BOOST_AUTO_TEST_CASE( cpt_mine_load_test )
{
    auto vertex1 = std::make_shared<bn::vertex_t>();
    auto vertex2 = std::make_shared<bn::vertex_t>();
    auto vertex3 = std::make_shared<bn::vertex_t>();

    vertex1->selectable_num = 2;
    vertex2->selectable_num = 2;
    vertex3->selectable_num = 3;

    bn::cpt_t cpt({vertex1, vertex2}, vertex3);

    bn::condition_t cond_a1 = { {vertex1, 0}, {vertex2, 0} };
    bn::condition_t cond_a2 = { {vertex1, 1}, {vertex2, 0} };
    bn::condition_t cond_a3 = { {vertex1, 0}, {vertex2, 1} };
    bn::condition_t cond_a4 = { {vertex1, 1}, {vertex2, 1} };

    cpt[cond_a1].second = {0.25,0.09,0.66};
    cpt[cond_a2].second = {0.10,0.15,0.75};
    cpt[cond_a3].second = {0.15,0.45,0.40};
    cpt[cond_a4].second = {0.00,0.50,0.50};

// OLD Style([[deprecated]])
//    bn::cpt_t cpt = {
//        {cond_a1, {0.25,0.09,0.66}}, {cond_a2, {0.10,0.15,0.75}},
//        {cond_a3, {0.15,0.45,0.40}}, {cond_a4, {0.00,0.50,0.50}}
//    };

    BOOST_CHECK(cpt[cond_a1].second == std::vector<double>({0.25,0.09,0.66}));
    BOOST_CHECK(cpt[cond_a2].second == std::vector<double>({0.10,0.15,0.75}));
    BOOST_CHECK(cpt[cond_a3].second == std::vector<double>({0.15,0.45,0.40}));
    BOOST_CHECK(cpt[cond_a4].second == std::vector<double>({0.00,0.50,0.50}));

    bn::condition_t filter1;
    bn::condition_t filter2 = {{vertex1, 0}};
    bn::condition_t filter3 = {{vertex2, 0}};
    bn::condition_t filter4 = {{vertex1, 1}, {vertex2, 1}};

    BOOST_CHECK(cpt.filter(filter1).size() == 4);
    BOOST_CHECK(cpt.filter(filter2).size() == 2);
    BOOST_CHECK(cpt.filter(filter3).size() == 2);
    BOOST_CHECK(cpt.filter(filter4).size() == 1);
}

BOOST_AUTO_TEST_CASE( cpt_other_load_test )
{
    auto vertex1 = std::make_shared<bn::vertex_t>();
    auto vertex2 = std::make_shared<bn::vertex_t>();
    auto vertex3 = std::make_shared<bn::vertex_t>();

    vertex1->selectable_num = 2;
    vertex2->selectable_num = 2;
    vertex3->selectable_num = 3;

    bn::cpt_t cpt({vertex1, vertex2}, vertex3);

    bn::condition_t cond_a1 = { {vertex1, 0}, {vertex2, 0} };
    bn::condition_t cond_a2 = { {vertex1, 1}, {vertex2, 0} };
    bn::condition_t cond_a3 = { {vertex1, 0}, {vertex2, 1} };
    bn::condition_t cond_a4 = { {vertex1, 1}, {vertex2, 1} };

    cpt[cond_a1].second = {0.25,0.09,0.66};
    cpt[cond_a2].second = {0.10,0.15,0.75};
    cpt[cond_a3].second = {0.15,0.45,0.40};
    cpt[cond_a4].second = {0.00,0.50,0.50};

// OLD Style([[deprecated]])
//    bn::cpt_t cpt = {
//        {cond_a1, {0.25,0.09,0.66}}, {cond_a2, {0.10,0.15,0.75}},
//        {cond_a3, {0.15,0.45,0.40}}, {cond_a4, {0.00,0.50,0.50}}
//    };

    bn::condition_t cond_b1 = { {vertex1, 0}, {vertex2, 1} };
    bn::condition_t cond_b2 = { {vertex1, 0}, {vertex2, 0} };
    bn::condition_t cond_b3 = { {vertex1, 1}, {vertex2, 0} };
    bn::condition_t cond_b4 = { {vertex1, 1}, {vertex2, 1} };
    BOOST_CHECK(cpt[cond_b1].second == std::vector<double>({0.15,0.45,0.40}));
    BOOST_CHECK(cpt[cond_b2].second == std::vector<double>({0.25,0.09,0.66}));
    BOOST_CHECK(cpt[cond_b3].second == std::vector<double>({0.10,0.15,0.75}));
    BOOST_CHECK(cpt[cond_b4].second == std::vector<double>({0.00,0.50,0.50}));
}

