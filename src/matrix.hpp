#ifndef BNI_MATRIX_HPP
#define BNI_MATRIX_HPP

#include <vector>
#include <boost/multi_array.hpp>

namespace bn {

// 制約:
// アクセスは y-x
using matrix_type = boost::multi_array<double, 2>;

// 行列同士・行列と整数の積(定義)
matrix_type operator*(matrix_type const& lhs, matrix_type const& rhs);
template<class T> matrix_type operator*(matrix_type const& rhs, T const& lhs);
template<class T> matrix_type operator*(T const& lhs, matrix_type const& rhs);

} // namespace bn

#include "matrix_impl.hpp"

#endif // #ifndef BNI_MATRIX_HPP

