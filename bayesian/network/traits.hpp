/**
* @file traits.hpp
* @brief Implementation of traits in bn::network.
* @author godai_0519
* @date 09/09/2016
*/

#ifndef BAYESIAN_NETWORKS_NETWORK_TRAITS_HPP
#define BAYESIAN_NETWORKS_NETWORK_TRAITS_HPP

#include <memory>
#include <type_traits>

namespace bn {
namespace traits {

//! Change shared_ptr<T> into shared_ptr<const T>.
template<class T>
struct add_const_shared {
    typedef std::shared_ptr<std::add_const_t<typename T::element_type>> type;
};

//! Change shared_ptr<T> into shared_ptr<const T>.
template<class T>
using add_const_shared_t = typename add_const_shared<T>::type;

} // namespace traits
} // namespace bn

#endif // BAYESIAN_NETWORKS_NETWORK_TRAITS_HPP
