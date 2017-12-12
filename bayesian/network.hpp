/**
* @file network.hpp
* @brief Provide a functions for each graph representing class.
* @author godai_0519
* @date 09/04/2016
*/

#ifndef BAYESIAN_NETWORKS_NETWORK_HPP
#define BAYESIAN_NETWORKS_NETWORK_HPP

#include <memory>
#include <vector>
#include <exception>
#include <bayesian/network/component.hpp>

namespace bn {

//! Provide a functions for each graph representing class.
/*! @tparam RepresentMethod: bn::adjacency_list or adjacency_matrix.
**/
template<class RepresentMethod>
class network {
public:
    using this_type = network<RepresentMethod>;
    using random_variable_ptr = component::random_variable_ptr;
    using node_ptr = component::node_ptr;
    using arc_ptr = component::arc_ptr;

    //! (Default ctor) Initialize nodes and arcs as empty.
    network()
        : pimpl_(new RepresentMethod())
    {
    }

    //! (Move ctor) Initialize nodes and arcs by moving a parameter.
    /*!        The parameter will be destroyed. */
    network(this_type&&) = default;

    //! (Default dtor)
    virtual ~network() = default;

    //! (Move operator=) Initialize nodes and arcs by copying a parameter.
    /*! The parameter will be destroyed. */
    this_type& operator=(this_type&&) = default;

    //! Clone this instance.
    /*! Cloned instance have the common pointers of nodes and arcs to this instance.

        @return      New cloned instance.
    **/
    this_type clone() const
    {
        std::unique_ptr<RepresentMethod> copied_ptr(new RepresentMethod(*pimpl_));
        return this_type(std::move(copied_ptr));
    }

    //! Register a node into network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @return      The node which is identical with the parameter.
    **/
    node_ptr add_node()
    {
        random_variable_ptr rv(new component::random_variable());
        node_ptr new_node(new component::node(rv));
        return pimpl_->add_node(new_node);
    }

    //! Register a cloned node of the argument into network.
    /*! An argument node and a return node are difference; however, they have the same instance of bn::random_variable.
        Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   node: a node to clone.
        @return      The node which is identical with the parameter.
    **/
    node_ptr add_clone_node(node_ptr const& node)
    {
        node_ptr ptr(new component::node(node->get()));
        return pimpl_->add_node(ptr);
    }

    //! Remove the argument node from network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   node: Target node.
        @return      true if the node was successfully removed; otherwise false.
    **/
    bool remove_node(node_ptr const& node)
    {
        return pimpl_->remove_node(node);
    }

    //! Remove all nodes from network.
    /*! Basic guarantee: if an exception is thrown, this instance is in a valid state. */
    void remove_all_node()
    {
        //TODO: there are some improvements
        for(auto const& node : pimpl_->all_node())
        {
            pimpl_->remove_node(node);
        }
    }

    //! Register an arc into network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   from: a node which is a start point of the arc.
        @param[in]   to: a node which is an end point of the arc.
        @return      The arc which is identical with the parameter.
    **/
    arc_ptr add_arc(node_ptr const& from, node_ptr const& to)
    {
        arc_ptr new_arc(new component::arc());
        return pimpl_->add_arc(new_arc, from, to);
    }

    //! Remove the argument arc from network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   arc: Target arc.
        @return      true if the node was successfully removed; otherwise false.
    **/
    bool remove_arc(arc_ptr const& arc)
    {
        return pimpl_->remove_arc(arc);
    }

    //! Remove the argument arc from network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   from: a node which is a start point of the arc wanted to remove.
        @param[in]   to: a node which is an end point of the arc wanted to remove.
        @return      true if the node was successfully removed; otherwise false.
    **/
    bool remove_arc(node_ptr const& from, node_ptr const& to)
    {
        return pimpl_->remove_arc(from, to);
    }

    //! Remove all arcs from network.
    /*! Basic guarantee: if an exception is thrown, this instance is in a valid state. */
    void remove_all_arc()
    {
        //TODO: there are some improvements
        for (auto const& arc : pimpl_->all_arc())
        {
            pimpl_->remove_arc(arc);
        }
    }

    //! Change a direction of argument arc.
    /*! Basic guarantee: if an exception is thrown, this instance is in a valid state.

        @return      true if the arc was successfully changed into an inverse one; otherwise false.
    **/
    bool change_direction(arc_ptr const& arc)
    {
        // Eraseable (that is, whether the arc exists or not)
        auto const new_source = pimpl_->target(arc);
        auto const new_target = pimpl_->source(arc);
        if(!new_source || !new_target) return false;

        // Try remove
        if(!pimpl_->remove_arc(arc)) return false;

        // Try addition of inversed arc
        if(!pimpl_->add_arc(arc, new_source, new_target))
        {
            // Restore
            if(!pimpl_->add_arc(arc, new_source, new_target))
            {
                //TODO: Libraries Exception
                throw std::runtime_error("");
            }
        }

        return true;
    }

    //! Check whether two nodes are adjacent by a arc in the network.
    /*! No-throw guarantee: never throws exceptions.

        @param[in]   from: a node which is a start point.
        @param[in]   to: a node which is an end point.
        @return      true if from and to are adjacent; otherwise false.
    **/
    arc_ptr is_adjacent(node_ptr const& from, node_ptr const& to) const noexcept
    {
        return pimpl_->is_adjacent(from, to);
    }

    //! Check whether a node and an arc are connected in the network.
    /*! No-throw guarantee: never throws exceptions.

        @param[in]   node: a node which is a start or end point.
        @param[in]   arc: an arc.
        @return      > 0 : the node is a start (parent) node of the arc;
                     < 0 : the node is an end (child) node of the arc;
                     = 0 : the node does not connect to the arc.
    **/
    int is_connect(node_ptr const& node, arc_ptr const& arc) const noexcept
    {
        return pimpl_->is_connect(node, arc);
    }

    //! Obtain a source (parent) node of an arc.
    /*! No-throw guarantee: never throws exceptions.

        @param[in]   arc: an arc
        @return      a source node pointer if it is exist; otherwise false
    **/
    node_ptr source(arc_ptr const& arc) const
    {
        return pimpl_->source(arc);
    }

    //! Obtain a target (child) node of an arc.
    /*! No-throw guarantee: never throws exceptions.

        @param[in]   arc: an arc.
        @return      a target node pointer if it is exist; otherwise false.
    **/
    node_ptr target(arc_ptr const& arc) const
    {
        return pimpl_->target(arc);
    }

    //! Obtain a list of all nodes which are parent of the child (argument).
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   child: a node which is an end point.
        @return      a list of all nodes.
    **/
    std::vector<node_ptr> parent_nodes(node_ptr const& node) const
    {
        return pimpl_->parent_nodes(node);
    }

    //! Obtain a list of all nodes which are child of the parent (argument).
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   parent: a node which is an start point.
        @return      a list of all nodes.
    **/
    std::vector<node_ptr> child_nodes(node_ptr const& node) const
    {
        return pimpl_->child_nodes(node);
    }

    //! Obtain all node pointers.
    std::vector<node_ptr> all_node() const
    {
        return pimpl_->all_node();
    }

    //! Obtain all arc pointers.
    std::vector<arc_ptr> all_arc() const
    {
        return pimpl_->all_arc();
    }

private:
    network(std::unique_ptr<RepresentMethod> pimpl)
        : pimpl_(std::move(pimpl))
    {
    }

    std::unique_ptr<RepresentMethod> pimpl_;
};

} // namespace bn

#endif // BAYESIAN_NETWORKS_NETWORK_HPP
