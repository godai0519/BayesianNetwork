/**
* @file component.hpp
* @brief aaaaa
* @author godai_0519
* @date 09/09/2016
*/

#ifndef BAYESIAN_NETWORKS_NETWORK_COMPONENT_HPP
#define BAYESIAN_NETWORKS_NETWORK_COMPONENT_HPP

#include <memory>

namespace bn {
namespace component {

//! @brief A class representing a random variable.
/*!        If max_value equals N, the random variable equals 0, 1,..., or N-1.

    @attention Default value of the max_value is equals to 0. You must set a value before using this.
**/
class random_variable {
public:
    //! The maximum value that this random variable can take.
    std::size_t max_value = 0;
};

//! @brief A class representing a node in network.
/*!        Each node has a corresponding random variable.
**/
class node {
public:
    using random_variable_ptr = std::shared_ptr<random_variable>;

    //! @brief Constructor unique to this class.
    /*! @param[in]   random_variable which is correnponding to this node. */
    explicit node(random_variable_ptr const& random_variable)
        : random_variable_(random_variable)
    {
    }

    //! @brief get the registered random variable.
    /*! @return registered random variable */
    random_variable_ptr get() const noexcept { return random_variable_; }

private:
    random_variable_ptr const random_variable_;
};


//! @brief A class representing a arc in network.
/*!        An arc cannot have any data. */
class arc {};

using random_variable_ptr = std::shared_ptr<random_variable>;
using node_ptr = std::shared_ptr<node>;
using arc_ptr = std::shared_ptr<arc>;

} // namespace component

} // namespace bn


#endif // BAYESIAN_NETWORKS_NETWORK_COMPONENT_HPP
