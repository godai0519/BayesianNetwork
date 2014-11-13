#ifndef BNI_MATRIX_HPP
#define BNI_MATRIX_HPP

#include <vector>
#include <boost/multi_array.hpp>

namespace bn {

// 制約:
// アクセスは y-x
typedef boost::multi_array<double, 2> matrix_type;


} // namespace bn

// 行列同士・行列と整数の積(定義)
bn::matrix_type operator*(bn::matrix_type const& lhs, bn::matrix_type const& rhs);
template<class T> bn::matrix_type operator*(bn::matrix_type const& rhs, T const& lhs);
template<class T> bn::matrix_type operator*(T const& lhs, bn::matrix_type const& rhs);

#include "matrix_impl.hpp"

#endif // #ifndef BNI_MATRIX_HPP

