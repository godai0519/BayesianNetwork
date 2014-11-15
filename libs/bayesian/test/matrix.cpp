#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "bayesian/matrix.hpp"

// teacher == mat1 * mat2 * mat3
BOOST_AUTO_TEST_CASE( operator_product )
{
    std::vector<double> const source1 = { 2, -1,  2, -1,  1,  0,  1,  0,  1};
    std::vector<double> const source2 = { 1,  0,  0,  0,  2,  0,  0,  0,  2};
    std::vector<double> const source3 = {-1, -1,  2, -1,  0,  2,  1,  1, -1};
    std::vector<double> const teacher = { 4,  2, -4, -1,  1,  2,  1,  1,  0};

    bn::matrix_type mat1(boost::extents[3][3]);
    bn::matrix_type mat2(boost::extents[3][3]);
    bn::matrix_type mat3(boost::extents[3][3]);
    bn::matrix_type matt(boost::extents[3][3]);

    mat1.assign(source1.begin(), source1.end());
    mat2.assign(source2.begin(), source2.end());
    mat3.assign(source3.begin(), source3.end());
    matt.assign(teacher.begin(), teacher.end());

    auto const result = mat1 * mat2 * mat3;
    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            // 誤差 0.00001%検査
            BOOST_CHECK_CLOSE(result[i][j], matt[i][j], 0.00001);
        }
    }
}

