/**
* @file network.hpp
* @brief aaaaa
* @author godai_0519
* @date 09/09/2016
*/

#ifndef BAYESIAN_NETWORKS_NETWORK_TRAITS_HPP
#define BAYESIAN_NETWORKS_NETWORK_TRAITS_HPP

#include <memory>
#include <type_traits>

namespace bn {
namespace traits {

template<class T>
struct add_const_shared {
    typedef std::shared_ptr<std::add_const_t<typename T::element_type>> type;
};

template<class T>
using add_const_shared_t = typename add_const_shared<T>::type;


} // namespace traits
} // namespace bn

#endif // BAYESIAN_NETWORKS_NETWORK_TRAITS_HPP
