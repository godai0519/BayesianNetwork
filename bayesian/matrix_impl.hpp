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
    for(matrix_type::size_type i = 0; i < result.shape()[0]; ++i)
    {
        for(matrix_type::size_type j = 0; j < result.shape()[1]; ++j)
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

