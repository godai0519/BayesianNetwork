#ifndef BNI_MATRIX_IMPL_HPP
#define BNI_MATRIX_IMPL_HPP

#include <utility>
#include "matrix.hpp"

namespace {
    using bn::matrix_type;
}

// 行列の整数倍
template<class T>
matrix_type operator*(matrix_type const& rhs, T const& lhs)
{
    // after copy, process all elem
    auto result = rhs;
    for(std::size_t i = 0; i < result.height(); ++i)
    {
        for(std::size_t j = 0; j < result.width(); ++j)
        {
            result[i][j] *= lhs;
        }
    }
    return result;
}

// delegate "integer * matrix" to "matrix * integer"
template<class T>
matrix_type operator*(T const& lhs, matrix_type const& rhs)
{
    return rhs * lhs;
}

#endif // #ifndef BNI_MATRIX_IMPL_HPP
