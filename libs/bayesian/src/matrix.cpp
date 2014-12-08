#include <vector>
#include "bayesian/matrix.hpp"

using bn::matrix_type;

matrix_type operator*(matrix_type const& lhs, matrix_type const& rhs)
{
    assert(lhs.width() == rhs.height()); // requirement for calculating the product

    matrix_type result(lhs.height(), rhs.width(), 0.0);
    for(std::size_t i = 0; i < result.height(); ++i)
    {
        for(std::size_t j = 0; j < result.width(); ++j)
        {
            for(std::size_t k = 0; k < lhs.width(); ++k)
            {
                result[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }

    return result;
}

