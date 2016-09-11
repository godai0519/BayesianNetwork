/**
* @file adjacency_matrix.hpp
* @brief Implementation of Adjacency Matrix using in network class
* @author godai_0519
* @date 09/09/2016
*/

#ifndef BAYESIAN_NETWORKS_NETWORK_ADJACENCY_MATRIX_HPP
#define BAYESIAN_NETWORKS_NETWORK_ADJACENCY_MATRIX_HPP

#include <algorithm> //find_if
#include <memory>
#include <unordered_map>
#include <vector>
#include "component.hpp"
#include "traits.hpp"

namespace bn {

//! Implement class of adjacency matrix for bn::network.
class adjacency_matrix {
public:
    using this_type = adjacency_matrix;
    using node_ptr = component::node_ptr;
    using arc_ptr = component::arc_ptr;
    using node_const_ptr = traits::add_const_shared_t<component::node_ptr>;
    using arc_const_ptr = traits::add_const_shared_t<component::arc_ptr>;

    using stored_node_type = std::list<node_ptr>;
    using node_dictionary_type = std::unordered_map<node_const_ptr, std::size_t>;
    using endpoint_dictionary_type =
        std::unordered_map<
        arc_const_ptr,
        std::pair<
        node_const_ptr,
        node_const_ptr>>;
    using adjacency_type = std::vector<std::vector<arc_ptr>>;

    //! (Default ctor) Initialize nodes and arcs as empty.
    adjacency_matrix() = default;

    //! (Copy ctor) Initialize nodes and arcs by copying a parameter.
    /*!        The parameter will "not" be destroyed (still valid). */
    adjacency_matrix(this_type const&) = default;

    //! (Move ctor) Initialize nodes and arcs by moving a parameter.
    /*!        The parameter will be destroyed. */
    adjacency_matrix(this_type&&) = default;

    //! @brief (Default dtor)
    virtual ~adjacency_matrix() = default;

    //! @brief (Copy operator=) Initialize nodes and arcs by copying a parameter.
    /*!        The parameter will "not" be destroyed (still valid). */
    this_type& operator=(this_type const& rhs)
    {
        *this = adjacency_matrix{ rhs };
        return *this;
    }

    //! @brief (Move operator=) Initialize nodes and arcs by copying a parameter.
    /*!        The parameter will be destroyed. */
    this_type& operator=(this_type&&) = default;

    //! Register a node into network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        param[in]   node: Target node.
        @return      The node which is identical with the parameter.
    **/
    node_ptr add_node(node_ptr const& node)
    {
        std::size_t const new_size = node_dic_.size() + 1;
        std::vector<arc_ptr> new_line(new_size);

        stored_node_.push_back(node);
        try
        {
            node_dic_[node] = new_size - 1;
            try
            {
                for(auto& vector : matrix_)
                    vector.resize(new_size);

                matrix_.push_back(std::move(new_line));
            }
            catch (std::exception const&)
            {
                for (auto& vector : matrix_)
                    vector.resize(new_size - 1); // No-throw guarantee

                node_dic_.erase(node);
                throw;
            }
        }
        catch(std::exception const&)
        {

            stored_node_.pop_back();
            throw;
        }

        return node;
    }


    //! Remove the argument node from network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   node: Target node.
        @return      true if the node was successfully removed; otherwise false.
    **/
    bool remove_node(node_ptr const& node)
    {
        auto it = std::find(stored_node_.crbegin(), stored_node_.crend(), node);
        if (it == stored_node_.crend()) return false; // Not registered

        node_dictionary_type backup_node_dic = node_dic_;
        adjacency_type backup_matrix = matrix_;

        auto const remove_index = node_dic_[node];
        try
        {
            matrix_.erase(matrix_.cbegin() + remove_index);
            for(auto& vector : matrix_)
                vector.erase(vector.cbegin() + remove_index);

            try
            {
                for(auto& elem : node_dic_)
                {
                    if(elem.second > remove_index) --elem.second;
                }
            }
            catch (std::exception const&)
            {
                node_dic_ = backup_node_dic;
                throw;
            }
        }
        catch(std::exception const&)
        {
            matrix_ = std::move(backup_matrix);
            throw;
        }

        for(auto it = endpoint_dic_.begin(); it != endpoint_dic_.end();)
        {
            if (it->second.first == node || it->second.second == node)
                it = endpoint_dic_.erase(it);
            else ++it;
        }

        stored_node_.erase((++it).base());
        return true;
    }

    //! Register an arc into network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   arc: an arc.
        @param[in]   from: a node which is a start point of the arc.
        @param[in]   to: a node which is an end point of the arc.
        @return      The arc which is identical with the parameter.
    **/
    arc_ptr add_arc(arc_ptr const& arc, node_ptr const& from, node_ptr const& to)
    {
        auto const from_index = node_dic_[from];
        auto const to_index = node_dic_[to];
        if (from_index >= matrix_.size() || to_index >= matrix_[from_index].size())
        {
            //TODO:
            throw std::runtime_error("");
        }

        auto endpoint = std::pair<node_const_ptr, node_const_ptr>(from, to);
        endpoint_dic_[arc] = std::move(endpoint);

        matrix_[from_index][to_index] = arc;
        return arc;
    }

    //! Remove the argument arc from network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   arc: Target arc.
        @return      true if the node was successfully removed; otherwise false.
    **/
    bool remove_arc(arc_ptr const& arc) noexcept
    {
        auto it = endpoint_dic_.find(arc);
        if (it == endpoint_dic_.cend()) return false; // Not registered

        return remove_arc(arc, it->second.first, it->second.second);
    }

    //! Remove the argument arc from network.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   from: a node which is a start point of the arc wanted to remove.
        @param[in]   to: a node which is an end point of the arc wanted to remove.
        @return      true if the node was successfully removed; otherwise false.
    **/
    bool remove_arc(node_ptr const& from, node_ptr const& to) noexcept
    {
        auto endpoint_it = std::find_if(
            endpoint_dic_.cbegin(), endpoint_dic_.cend(),
            [&from, &to](auto const& elem) {
            return elem.second.first == from && elem.second.second == to;
        });
        if (endpoint_it == endpoint_dic_.cend()) return false; // Not registered

        return remove_arc(std::const_pointer_cast<component::arc>(endpoint_it->first), from, to);
    }

    //! Check whether two nodes are adjacent by a arc in the network.
    /*! No-throw guarantee: never throws exceptions.

        @param[in]   from: a node which is a start point.
        @param[in]   to: a node which is an end point.
        @return      true if from and to are adjacent; otherwise false.
    **/
    arc_ptr is_adjacent(node_ptr const& from, node_ptr const& to) const
    {
        auto const from_index = node_dic_.at(from);
        auto const to_index = node_dic_.at(to);

        if(from_index >= matrix_.size() || to_index >= matrix_[from_index].size()) return nullptr;
        return matrix_[from_index][to_index];
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
        auto endpoint_it = endpoint_dic_.find(arc);
        if(endpoint_it == endpoint_dic_.cend()) return 0; // Not registered


        return endpoint_it->second.first  == node ?  1 :
               endpoint_it->second.second == node ? -1 :
                                                     0 ;
    }

    //! Obtain a source (parent) node of an arc.
    /*! No-throw guarantee: never throws exceptions.

        @param[in]   arc: an arc
        @return      a source node pointer if it is exist; otherwise false
    **/
    node_ptr source(arc_ptr const& arc) const noexcept
    {
        auto endpoint_it = endpoint_dic_.find(arc);
        if (endpoint_it == endpoint_dic_.cend()) return nullptr; // Not registered

        return std::const_pointer_cast<component::node>(endpoint_it->second.first);
    }

    //! Obtain a target (child) node of an arc.
    /*! No-throw guarantee: never throws exceptions.

        @param[in]   arc: an arc.
        @return      a target node pointer if it is exist; otherwise false.
    **/
    node_ptr target(arc_ptr const& arc) const noexcept
    {
        auto endpoint_it = endpoint_dic_.find(arc);
        if (endpoint_it == endpoint_dic_.cend()) return nullptr; // Not registered

        return std::const_pointer_cast<component::node>(endpoint_it->second.second);
    }

    //! Obtain a list of all nodes which are parent of the child (argument).
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   child: a node which is an end point.
        @return      a list of all nodes.
    **/
    std::vector<node_ptr> parent_nodes(node_ptr const& child) const
    {
        auto const child_index = node_dic_.at(child);

        std::vector<node_ptr> parents;
        parents.reserve(stored_node_.size());

        for(auto it = stored_node_.cbegin(); it != stored_node_.cend(); ++it)
        {
            auto const it_index = node_dic_.at(*it);
            if (matrix_[it_index][child_index])
                parents.push_back(*it);
        }

        parents.shrink_to_fit();
        return parents;
    }

    //! Obtain a list of all nodes which are child of the parent (argument).
    /*! Strong guarantee: if an exception is thrown, there are no changes in the network.

        @param[in]   parent: a node which is an start point.
        @return      a list of all nodes.
    **/
    std::vector<node_ptr> child_nodes(node_ptr const& parent) const
    {
        auto const parent_index = node_dic_.at(parent);

        std::vector<node_ptr> children;
        children.reserve(stored_node_.size());

        for(auto it = stored_node_.cbegin(); it != stored_node_.cend(); ++it)
        {
            auto const it_index = node_dic_.at(*it);
            if(matrix_[parent_index][it_index])
                children.push_back(*it);
        }

        children.shrink_to_fit();
        return children;
    }

    //! Obtain all node pointers.
    std::vector<node_ptr> all_node() const
    {
        return std::vector<node_ptr>(stored_node_.cbegin(), stored_node_.cend());
    }

    //! Obtain all arc pointers.
    std::vector<arc_ptr> all_arc() const
    {
        std::vector<arc_ptr> arcs;
        arcs.reserve(endpoint_dic_.size());

        for(auto const& endpoint : endpoint_dic_)
            arcs.push_back(std::const_pointer_cast<component::arc>(endpoint.first));

        arcs.shrink_to_fit();
        return arcs;
    }

private:
    bool remove_arc(arc_ptr const& arc, node_const_ptr const& from, node_const_ptr const& to) noexcept
    {
        auto const from_index = node_dic_[from];
        auto const to_index = node_dic_[to];
        if(from_index >= matrix_.size() || to_index >= matrix_[from_index].size()) return false;

        matrix_[from_index][to_index] = nullptr;
        endpoint_dic_.erase(arc);
        return true;
    }

    stored_node_type stored_node_;
    node_dictionary_type node_dic_;
    endpoint_dictionary_type endpoint_dic_;
    adjacency_type matrix_;
};

} // namespace bn

#endif // BAYESIAN_NETWORKS_NETWORK_ADJACENCY_MATRIX_HPP
