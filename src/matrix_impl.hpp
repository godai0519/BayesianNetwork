#ifndef BNI_MATRIX_IMPL_HPP
#define BNI_MATRIX_IMPL_HPP

#include <utility>
#include "matrix.hpp"

namespace bn {

// s—ñ‚Ì®””{
template<class T>
matrix_type operator*(matrix_type const& rhs, T const& lhs)
{
    // after copy, process all elem
    auto result = rhs;
    for(auto it = result.begin(); it != result.end(); ++it)
    {
        *it *= lhs;
    }
    return result;
}

// delegate "integer * matrix" to "matrix * integer" 
template<class T>
matrix_type operator*(T const& lhs, matrix_type const& rhs)
{
    return rhs * lhs;
}


} // namespace bn

#endif // #ifndef BNI_MATRIX_IMPL_HPP

