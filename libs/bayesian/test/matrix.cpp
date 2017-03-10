#define BOOST_TEST_MAIN
#include <array>
#include <random>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include "bayesian/matrix.hpp"

BOOST_AUTO_TEST_CASE(assign_test)
{
    // Generate m19937
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);
    std::generate(v.begin(), v.end(), std::ref(rnd));
    std::seed_seq seq(v.begin(), v.end());
    std::mt19937 engine(seq);
    std::uniform_int_distribution<std::size_t> int_dist(1, 200);
    std::uniform_real_distribution<double> real_dist(0.0, 1000);

    // Determine size
    std::size_t data_size = 1;
    std::vector<std::size_t> size(3);
    for(std::size_t i = 0; i < 3; ++i)
    {
        size[i] = int_dist(engine);
        data_size *= size[i];
    }

    // Containers
    bn::matrix<double> mat1, mat2, mat3;
    std::vector<double> data(data_size);
    std::generate(data.begin(), data.end(), [&real_dist, &engine]() { return real_dist(engine); });

    // Assign
    double const const_data = 100;
    mat1.assign(size, const_data);
    mat2.assign(size, data);
    mat3.assign(size, data.begin(), data.end());

    for(std::size_t i = 0; i < size[0]; ++i)
    {
        for(std::size_t j = 0; j < size[1]; ++j)
        {
            for (std::size_t k = 0; k < size[2]; ++k)
            {
                auto const value1 = mat1[std::array<std::size_t, 3>{{i, j, k}}];
                auto const value2 = mat2[std::array<std::size_t, 3>{{i, j, k}}];
                auto const value3 = mat3[std::array<std::size_t, 3>{{i, j, k}}];

                BOOST_CHECK_CLOSE(value1, const_data, 0.0001);
                BOOST_CHECK_CLOSE(value2, data[(i * size[1] + j) * size[2] + k], 0.0001);
                BOOST_CHECK_CLOSE(value3, data[(i * size[1] + j) * size[2] + k], 0.0001);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(resize)
{
    // Generate m19937
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);
    std::generate(v.begin(), v.end(), std::ref(rnd));
    std::seed_seq seq(v.begin(), v.end());
    std::mt19937 engine(seq);
    std::uniform_int_distribution<std::size_t> int_dist(5, 200);
    std::uniform_real_distribution<double> real_dist(0.0, 1000);

    // Determine size
    std::size_t data_size = 1;
    std::vector<std::size_t> size(3);
    for (std::size_t i = 0; i < 3; ++i)
    {
        size[i] = int_dist(engine);
        data_size *= size[i];
    }

    // Determine data
    std::vector<double> data(data_size);
    std::generate(data.begin(), data.end(), [&real_dist, &engine]() { return real_dist(engine); });

    bn::matrix<double> mat_original(size, data),
                       mat1(size, data), // Increase Dimension
                       mat2(size, data), // Decrease Dimension
                       mat3(size, data), // Increase Size
                       mat4(size, data); // Decrease Size

    // Resize
    mat1.resize(std::vector<std::size_t>{size[0], size[1], size[2], 5});
    mat2.resize(std::vector<std::size_t>{size[0], size[1]});
    mat3.resize(std::vector<std::size_t>{size[0] + 3, size[1] + 3, size[2] + 3});
    mat4.resize(std::vector<std::size_t>{size[0] - 3, size[1] - 3, size[2] - 3});

    for(std::size_t i = 0; i < size[0]; ++i)
    {
        for(std::size_t j = 0; j < size[1]; ++j)
        {
            double const elem_mat2 = mat2[std::array<std::size_t, 2>{{i, j}}];
            double const elem_cutz = mat_original[std::array<std::size_t, 3>{{i, j, 0}}];
            BOOST_CHECK_CLOSE(elem_mat2, elem_cutz, 0.0001);

            for(std::size_t k = 0; k < size[2]; ++k)
            {
                double const elem_original = mat_original[std::array<std::size_t, 3>{ {i, j, k}}];

                double const elem_mat1 = mat1[std::array<std::size_t, 4>{{i, j, k, 0}}];
                double const elem_mat3 = mat3[std::array<std::size_t, 3>{{i, j, k}}];
                BOOST_CHECK_CLOSE(elem_mat1, elem_original, 0.0001);
                BOOST_CHECK_CLOSE(elem_mat3, elem_original, 0.0001);

                if(i < size[0] - 3 && j < size[1] - 3 && k < size[2] - 3)
                {
                    double const elem_mat4 = mat4[std::array<std::size_t, 3>{ {i, j, k}}];
                    BOOST_CHECK_CLOSE(elem_mat4, elem_original, 0.0001);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(resize_default_value)
{
    // Generate m19937
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);
    std::generate(v.begin(), v.end(), std::ref(rnd));
    std::seed_seq seq(v.begin(), v.end());
    std::mt19937 engine(seq);
    std::uniform_int_distribution<std::size_t> int_dist(5, 200);
    std::uniform_real_distribution<double> real_dist(0.0, 1000);

    // Determine size
    std::size_t data_size = 1;
    std::vector<std::size_t> size(3);
    for (std::size_t i = 0; i < 3; ++i)
    {
        size[i] = int_dist(engine);
        data_size *= size[i];
    }

    // Determine data
    std::vector<double> data(data_size);
    std::generate(data.begin(), data.end(), [&real_dist, &engine]() { return real_dist(engine); });

    bn::matrix<double> mat_original(size, data),
        mat1(size, data), // Increase Dimension
        mat2(size, data), // Decrease Dimension
        mat3(size, data), // Increase Size
        mat4(size, data); // Decrease Size

    // Resize
    double const default_value = real_dist(engine);
    mat1.resize(std::vector<std::size_t>{size[0], size[1], size[2], 5}, default_value);
    mat2.resize(std::vector<std::size_t>{size[0], size[1]}, default_value);
    mat3.resize(std::vector<std::size_t>{size[0] + 3, size[1] + 3, size[2] + 3}, default_value);
    mat4.resize(std::vector<std::size_t>{size[0] - 3, size[1] - 3, size[2] - 3}, default_value);

    for (std::size_t i = 0; i < size[0]; ++i)
    {
        for (std::size_t j = 0; j < size[1]; ++j)
        {
            double const elem_mat2 = mat2[std::array<std::size_t, 2>{ {i, j}}];
            double const elem_cutz = mat_original[std::array<std::size_t, 3>{ {i, j, 0}}];
            BOOST_CHECK_CLOSE(elem_mat2, elem_cutz, 0.0001);

            for (std::size_t k = 0; k < size[2]; ++k)
            {
                double const elem_original = mat_original[std::array<std::size_t, 3>{ {i, j, k}}];

                double const elem_mat1 = mat1[std::array<std::size_t, 4>{ {i, j, k, 0}}];
                double const elem_mat3 = mat3[std::array<std::size_t, 3>{ {i, j, k}}];
                BOOST_CHECK_CLOSE(elem_mat1, elem_original, 0.0001);
                BOOST_CHECK_CLOSE(elem_mat3, elem_original, 0.0001);

                if (i < size[0] - 3 && j < size[1] - 3 && k < size[2] - 3)
                {
                    double const elem_mat4 = mat4[std::array<std::size_t, 3>{ {i, j, k}}];
                    BOOST_CHECK_CLOSE(elem_mat4, elem_original, 0.0001);
                }
            }
        }
    }


    for(std::size_t i = 0; i < size[0]; ++i)
    {
        for(std::size_t j = 0; j < size[1]; ++j)
        {
            for(std::size_t k = 0; k < size[2]; ++k)
            {
                for(std::size_t l = 1; l < 5; ++l)
                {
                    double const elem_mat1 = mat1[std::array<std::size_t, 4>{ {i, j, k, l}}];
                    BOOST_CHECK_CLOSE(elem_mat1, default_value, 0.0001);
                }
            }
        }
    }

    for(std::size_t i = 0; i < size[0] + 3; ++i)
    {
        for(std::size_t j = 0; j < size[1] + 3; ++j)
        {
            for(std::size_t k = 0; k < size[2] + 3; ++k)
            {
                if(i >= size[0] || j >= size[1] || k >= size[2])
                {
                    double const elem_mat3 = mat3[std::array<std::size_t, 4>{ {i, j, k}}];
                    BOOST_CHECK_CLOSE(elem_mat3, default_value, 0.0001);

                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(shrink_to_fit)
{
    // Generate m19937
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);
    std::generate(v.begin(), v.end(), std::ref(rnd));
    std::seed_seq seq(v.begin(), v.end());
    std::mt19937 engine(seq);
    std::uniform_int_distribution<std::size_t> int_dist(5, 200);
    std::uniform_real_distribution<double> real_dist(0.0, 1000);

    // Determine size
    std::size_t data_size = 1;
    std::vector<std::size_t> size(4);
    for (std::size_t i = 0; i < 4; ++i)
    {
        size[i] = int_dist(engine);
        data_size *= size[i];
    }

    // Determine data
    std::vector<double> data(data_size);
    std::generate(data.begin(), data.end(), [&real_dist, &engine]() { return real_dist(engine); });

    // Original Matrix
    bn::matrix<double> mat_original(size, data);
    mat_original.resize(std::vector<std::size_t>{size[0], 1, size[2], 1});

    // Shrink
    bn::matrix<double> mat_shrink = mat_original;
    mat_shrink.shrink_to_fit();

    BOOST_CHECK(mat_original.dims() == mat_shrink.dims());
    BOOST_CHECK(mat_original.sizes() == mat_shrink.sizes());
    BOOST_CHECK(mat_original.data().size() / data[1] / data[3] == mat_shrink.data().size());

    for(std::size_t i = 0; i < size[0]; ++i)
    {
        for(std::size_t j = 0; j < size[2]; ++j)
        {
            double value_original = mat_original[std::array<std::size_t, 4>{ {i, 1, j, 1}}];
            double value_shrink = mat_shrink[std::array<std::size_t, 4>{{i, 1, j, 1}}];
            BOOST_CHECK_CLOSE(value_original, value_shrink, 0.0001);
        }
    }
}

BOOST_AUTO_TEST_CASE(copy_move_test)
{
    // Generate m19937
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);
    std::generate(v.begin(), v.end(), std::ref(rnd));
    std::seed_seq seq(v.begin(), v.end());
    std::mt19937 engine(seq);
    std::uniform_int_distribution<std::size_t> int_dist(1, 200);
    std::uniform_real_distribution<double> real_dist(0.0, 1000);

    // Determine size & data
    std::size_t data_size = 1;
    std::vector<std::size_t> size(3);
    for (std::size_t i = 0; i < 3; ++i)
    {
        size[i] = int_dist(engine);
        data_size *= size[i];
    }
    std::vector<double> data(data_size);
    std::generate(data.begin(), data.end(), [&real_dist, &engine]() { return real_dist(engine); });

    // Copy
    bn::matrix<double> mat(size, data);
    bn::matrix<double> mat_copy = mat;
    BOOST_CHECK(mat.dims() == mat_copy.dims());
    BOOST_CHECK(mat.sizes() == mat_copy.sizes());
    BOOST_CHECK(mat.capacities() == mat_copy.capacities());
    BOOST_CHECK(mat.steps() == mat_copy.steps());
    BOOST_CHECK(mat.data() == mat_copy.data());

    bn::matrix<double> mat_move = std::move(mat);
    BOOST_CHECK(mat_move.dims() == mat_copy.dims());
    BOOST_CHECK(mat_move.sizes() == mat_copy.sizes());
    BOOST_CHECK(mat_move.capacities() == mat_copy.capacities());
    BOOST_CHECK(mat_move.steps() == mat_copy.steps());
    BOOST_CHECK(mat_move.data() == mat_copy.data());
}

BOOST_AUTO_TEST_CASE(dottest)
{
    bn::matrix<double> mat1(std::vector<std::size_t>{4}, std::vector<double>{1, 2, 3, 4});
    bn::matrix<double> mat2(std::vector<std::size_t>{4}, std::vector<double>{10, 20, 30, 40});
    bn::matrix<double> mat3(std::vector<std::size_t>{5}, std::vector<double>{10, 20, 30, 40, 50});
    bn::matrix<double> mat4(
        std::vector<std::size_t>{4, 4},
        std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
    );

    BOOST_CHECK(mat1.dot(mat2) == 300);
    BOOST_CHECK(mat1.dot(mat2) == mat2.dot(mat1));
    BOOST_CHECK(mat1.dot(mat2) == dot(mat1, mat2));

    BOOST_CHECK_THROW(dot(mat1, mat3), std::runtime_error);
    BOOST_CHECK_THROW(dot(mat1, mat4), std::runtime_error);
}

BOOST_DATA_TEST_CASE(multiplication_mm, boost::unit_test::data::xrange(10))
{
    // Generate m19937
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);
    std::generate(v.begin(), v.end(), std::ref(rnd));
    std::seed_seq seq(v.begin(), v.end());
    std::mt19937 engine(seq);
    std::uniform_int_distribution<std::size_t> int_dist(1, 200);
    std::uniform_real_distribution<double> real_dist(0.0, 1000);

    // Generate parameter
    std::size_t const row1 = int_dist(engine);
    std::size_t const col1 = int_dist(engine);
    std::size_t const col2 = int_dist(engine);

    // Initialize matrix
    bn::matrix<double> mat1(std::vector<std::size_t>{row1, col1});
    bn::matrix<double> mat2(std::vector<std::size_t>{col1, col2});
    for(auto& elem : mat1.data()) const_cast<double&>(elem) = real_dist(engine);
    for(auto& elem : mat2.data()) const_cast<double&>(elem) = real_dist(engine);

    auto const result = mat1 * mat2;
    for(std::size_t i = 0; i < row1; ++i)
    {
        for(std::size_t j = 0; j < col2; ++j)
        {
            double value = 0;
            for(std::size_t k = 0; k < col1; ++k)
                value += mat1[std::array<std::size_t, 2>{{i, k}}] * mat2[std::array<std::size_t, 2>{{k, j}}];

            double const calculated = result[std::array<std::size_t, 2>{ {i, j}}];
            BOOST_CHECK_CLOSE(calculated, value, 0.0001);
        }
    }
}

BOOST_DATA_TEST_CASE(multiplication_ms, boost::unit_test::data::xrange(10))
{
    // Generate m19937
    std::random_device rnd;
    std::vector<std::uint_least32_t> v(10);
    std::generate(v.begin(), v.end(), std::ref(rnd));
    std::seed_seq seq(v.begin(), v.end());
    std::mt19937 engine(seq);
    std::uniform_int_distribution<std::size_t> int_dist(1, 200);
    std::uniform_real_distribution<double> real_dist(0.0, 1000);

    // Generate parameter
    std::size_t const row = int_dist(engine);
    std::size_t const col = int_dist(engine);

    // Initialize matrix
    double scala = real_dist(engine);
    bn::matrix<double> mat(std::vector<std::size_t>{row, col});
    for(auto& elem : mat.data()) const_cast<double&>(elem) = real_dist(engine);

    auto const result1 = scala * mat;
    auto const result2 = mat * scala;
    for(std::size_t i = 0; i < row; ++i)
    {
        for(std::size_t j = 0; j < col; ++j)
        {
            double const original = mat[std::array<std::size_t, 2>{ {i, j}}];
            double const calculated1 = result1[std::array<std::size_t, 2>{{i, j}}];
            double const calculated2 = result2[std::array<std::size_t, 2>{{i, j}}];
            BOOST_CHECK_CLOSE(calculated1, original * scala, 0.0001);
            BOOST_CHECK_CLOSE(calculated1, calculated2, 0.0001);
        }
    }
}
