#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include "bayesian/cpt.hpp"

bn::cpt::rv_ptr generate_target(
	std::mt19937& mt,
	std::uniform_int_distribution<int>& max_value_dist
)
{
	bn::cpt::rv_ptr target_rv(new bn::component::random_variable());
	target_rv->max_value = max_value_dist(mt);
	return target_rv;
}

std::vector<bn::cpt::rv_ptr> generate_parents(
	std::mt19937& mt,
	std::uniform_int_distribution<int>& max_value_dist,
	std::size_t const parent_num
)
{
	std::vector<bn::cpt::rv_ptr> parent_rv(parent_num);
	std::generate(parent_rv.begin(), parent_rv.end(),
		[&mt, &max_value_dist]() -> auto {
		bn::cpt::rv_ptr rv(new bn::component::random_variable());
		rv->max_value = max_value_dist(mt);
		return rv;
	});

	return parent_rv;
}

bn::cpt::condition_type generate_condition(
	std::mt19937& mt,
	std::vector<bn::cpt::rv_ptr> parent_rv
)
{
	bn::cpt::condition_type cond;
	for (auto const& rv : parent_rv)
	{
		std::uniform_int_distribution<int> parent_dist(0, rv->max_value - 1);
		cond[rv] = parent_dist(mt);
	}
	return cond;
}

BOOST_DATA_TEST_CASE(cpt_zero_dimension, boost::unit_test::data::xrange(50), i)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> max_value_dist(1, 10);
	std::uniform_real_distribution<double> probability_dist(0, 1);

	auto const target_rv = generate_target(mt, max_value_dist);
	
    bn::cpt cpt(target_rv);
	bn::cpt::condition_type cond = {};
	bn::cpt::condition_type cond_another = { {target_rv, 1} };
	BOOST_CHECK(cpt.is_valid(cond));
	BOOST_CHECK(cpt.is_valid(cond_another));

    auto const indexes = cpt.condition_to_index(cond);
    BOOST_CHECK(indexes.size() == 1);
    BOOST_CHECK(indexes.at(0) == 0);

    auto& elem1 = cpt.at(cond);
    auto& elem2 = cpt[cond];
    BOOST_CHECK(elem1.size() == target_rv->max_value);
    BOOST_CHECK(elem2.size() == target_rv->max_value);

    std::vector<double> sample(target_rv->max_value);
    for(std::size_t j = 0; j < target_rv->max_value; ++j)
        sample[j] = probability_dist(mt);
    std::copy(sample.begin(), sample.end(), elem1.begin());
    for(std::size_t j = 0; j < target_rv->max_value; ++j)
        BOOST_CHECK_CLOSE(elem1[j], elem2[j], 0.0001);
}

BOOST_DATA_TEST_CASE(cpt_some_dimensions, boost::unit_test::data::xrange(50), parent_num)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> max_value_dist(1, 10);
	std::uniform_real_distribution<double> probability_dist(0, 1);

	auto const target_rv = generate_target(mt, max_value_dist);
	auto const parent_rv = generate_parents(mt, max_value_dist, ((parent_num + 1) % 10) + 1);
	
	bn::cpt cpt(target_rv, parent_rv);
	auto const cond = generate_condition(mt, parent_rv);
	auto cond_invalid = cond; cond_invalid.erase(cond_invalid.begin());
	BOOST_CHECK(cpt.is_valid(cond));
	BOOST_CHECK(cpt.is_valid(cond_invalid) == false);

	auto const indexes = cpt.condition_to_index(cond);
	BOOST_CHECK(indexes.size() == parent_rv.size());
	for(int j = 0; j < parent_rv.size(); ++j)
		BOOST_CHECK(indexes.at(j) == cond.at(parent_rv[j]));

	auto& elem1 = cpt.at(cond);
	auto& elem2 = cpt[cond];
	BOOST_CHECK(elem1.size() == target_rv->max_value);
	BOOST_CHECK(elem2.size() == target_rv->max_value);

	std::vector<double> sample(target_rv->max_value);
	for (std::size_t j = 0; j < target_rv->max_value; ++j)
		sample[j] = probability_dist(mt);
	std::copy(sample.begin(), sample.end(), elem1.begin());
	for (std::size_t j = 0; j < target_rv->max_value; ++j)
		BOOST_CHECK_CLOSE(elem1[j], elem2[j], 0.0001);
}

BOOST_DATA_TEST_CASE(cpt_copy, boost::unit_test::data::xrange(50), parent_num)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> max_value_dist(1, 10);
	std::uniform_real_distribution<double> probability_dist(0, 1);

	auto const target_rv = generate_target(mt, max_value_dist);
	auto const parent_rv = generate_parents(mt, max_value_dist, ((parent_num + 1) % 10) + 1);

	bn::cpt cpt(target_rv, parent_rv);
	auto const cond = generate_condition(mt, parent_rv);
	BOOST_CHECK(cpt.is_valid(cond));

	// Setting data
	auto& elem = cpt[cond];
	for (std::size_t j = 0; j < target_rv->max_value; ++j)
		elem[j] = probability_dist(mt);

	// Copy!
	bn::cpt cpt_copied = cpt;
	BOOST_CHECK(cpt_copied.is_valid(cond));

	// Comparing the target node
	BOOST_CHECK(cpt.rv() == cpt_copied.rv());

	// Comparing the parent nodes
	for(int j = 0; j < cpt.parents().size(); ++j)
		BOOST_CHECK(cpt.parents()[j] == cpt_copied.parents()[j]);

	// Comparing the probabilities of target node
	auto& elem_copied = cpt_copied[cond];
	for (std::size_t j = 0; j < target_rv->max_value; ++j)
		BOOST_CHECK_CLOSE(elem[j], elem_copied[j], 0.0001);
}

BOOST_AUTO_TEST_CASE(cpt_manager_for_same_rv)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> parent_num_dist(1, 10);
	std::uniform_int_distribution<int> max_value_dist(1, 10);
	std::uniform_real_distribution<double> probability_dist(0, 1);

	auto const target_rv1 = generate_target(mt, max_value_dist);
	auto const parent_rv1 = generate_parents(mt, max_value_dist, parent_num_dist(mt));
	bn::cpt cpt1(target_rv1, parent_rv1);
	bn::cpt cpt2(target_rv1, parent_rv1);

	auto const target_rv2 = generate_target(mt, max_value_dist);
	auto const parent_rv2 = generate_parents(mt, max_value_dist, parent_num_dist(mt));
	bn::cpt cpt3(target_rv2, parent_rv2);

	// Different nodes for the same r.v. and another node.
	auto const node1 = std::make_shared<bn::component::node>(target_rv1);
	auto const node2 = std::make_shared<bn::component::node>(target_rv1);
	auto const node3 = std::make_shared<bn::component::node>(target_rv2);

	bn::cpt_manager manager;
	manager.enroll(node1, cpt1);
	manager.enroll(node2, std::move(cpt2));
	manager.enroll(node3, cpt3);

	BOOST_CHECK(&manager.at(node1) != &manager[node2]);
	BOOST_CHECK(&manager.at(node1) != &manager[node3]);

	manager.unenroll(node1);
	manager.unenroll(node2);
	manager.unenroll(node3);

	BOOST_CHECK_THROW(manager.at(node1), std::out_of_range);
	BOOST_CHECK_THROW(manager.at(node2), std::out_of_range);
	BOOST_CHECK_THROW(manager.at(node3), std::out_of_range);
}
