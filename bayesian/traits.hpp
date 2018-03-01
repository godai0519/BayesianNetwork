/**
* @file traits.hpp
* @brief Implementation of traits in bn::network.
* @author godai_0519
* @date 09/09/2016
*/

#ifndef BAYESIAN_NETWORKS_TRAITS_HPP
#define BAYESIAN_NETWORKS_TRAITS_HPP

#include <iterator>
#include <type_traits>
#include <utility>

namespace bn {
namespace traits {

#if __cpp_lib_void_t
    using std::void_t;
#else
    template<class...>  using void_t = void;
#endif

namespace is_iterable_impl { // {{{

using std::begin;
using std::end;

template<class, class = void>
struct enable_begin_end : std::false_type {};

template<class T>
struct enable_begin_end<T,
    void_t<decltype(begin(std::declval<const T&>()), end(std::declval<const T&>()))>>
    : std::true_type
{};

template < typename T >
constexpr bool enable_begin_end_v = enable_begin_end<T>::value;

} // namespace is_iterable_impl }}}

template < typename T >
struct is_range : std::integral_constant<bool,
    is_iterable_impl::enable_begin_end_v<T> ||
    std::is_array<T>::value>
{};

template < typename T >
constexpr bool is_range_v = is_range<T>::value;

} // namespace traits
} // namespace bn

#endif // BAYESIAN_NETWORKS_TRAITS_HPP
