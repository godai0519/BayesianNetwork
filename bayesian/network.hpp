/**
* @file network.hpp
* @brief aaaaa
* @author godai_0519
* @date 09/04/2016
*/

#ifndef BAYESIAN_NETWORKS_NETWORK_HPP
#define BAYESIAN_NETWORKS_NETWORK_HPP

#include <memory>
#include <vector>
#include <exception>
#include "network/component.hpp"

namespace bn {

template<class RepresentMethod>
class network {
public:
    using this_type = network<RepresentMethod>;
    using random_variable_ptr = component::random_variable_ptr;
    using node_ptr = component::node_ptr;
    using arc_ptr = component::arc_ptr;

    // default ctor
    // Initialize nodes and arcs as empty
    network()
        : pimpl_(new RepresentMethod())
    {
    }

    // move ctor
    network(this_type&&) = default;

    // default dtor
    virtual ~network() = default;

    // move operator=
    this_type& operator=(this_type&&) = default;

    this_type clone()
    {
        std::unique_ptr<RepresentMethod> copied_ptr(new RepresentMethod(*pimpl_));
        return this_type(std::move(copied_ptr));
    }

    node_ptr add_node()
    {
        random_variable_ptr rv(new component::random_variable());
        node_ptr new_node(new component::node(rv));
        return pimpl_->add_node(new_node);
    }

    node_ptr add_clone_node(node_ptr const& node)
    {
        node_ptr ptr(new component::node(node->get()));
        return pimpl_->add_node(ptr);
    }

    bool remove_node(node_ptr const& node)
    {
        return pimpl_->remove_node(node);
    }

    void remove_all_node()
    {
        //TODO: there are some improvements
        for(auto const& node : pimpl_->all_node())
        {
            pimpl_->remove_node(node);
        }
    }

    arc_ptr add_arc(node_ptr const& from, node_ptr const& to)
    {
        arc_ptr new_arc(new component::arc());
        return pimpl_->add_arc(new_arc, from, to);
    }

    bool remove_arc(node_ptr const& from, node_ptr const& to)
    {
        return pimpl_->remove_arc(from, to);
    }

    bool remove_arc(arc_ptr const& arc)
    {
        return pimpl_->remove_arc(arc);
    }

    void remove_all_arc()
    {
        //TODO: there are some improvements
        for (auto const& arc : pimpl_->all_arc())
        {
            pimpl_->remove_arc(arc);
        }
    }

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

    arc_ptr is_adjacent(node_ptr const& from, node_ptr const& to) const noexcept
    {
        return pimpl_->is_adjacent(from, to);
    }

    int is_connect(node_ptr const& node, arc_ptr const& arc) const noexcept
    {
        return pimpl_->is_connect(node, arc);
    }

    node_ptr source(arc_ptr const& arc) const
    {
        return pimpl_->source(arc);
    }

    node_ptr target(arc_ptr const& arc) const
    {
        return pimpl_->target(arc);
    }

    std::vector<node_ptr> parent_nodes(node_ptr const& node) const
    {
        return pimpl_->parent_nodes(node);
    }

    std::vector<node_ptr> child_nodes(node_ptr const& node) const
    {
        return pimpl_->child_nodes(node);
    }

    std::vector<node_ptr> all_node() const
    {
        return pimpl_->all_node();
    }

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
