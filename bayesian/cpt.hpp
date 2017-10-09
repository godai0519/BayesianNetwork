/**
* @file cpt.hpp
* @brief A wrapper class of bn::matrix provides access methods to the CPT (Conditional Probability Table).
* @author godai_0519
* @date 12/01/2016
*/

#ifndef BAYESIAN_NETWORKS_CPT_HPP
#define BAYESIAN_NETWORKS_CPT_HPP

#include <unordered_map>
#include <vector>
#include "matrix.hpp"
#include "network/component.hpp"

namespace bn {

//! A wrapper class of bn::matrix provides access methods to the CPT (Conditional Probability Table).
/*! The CPT is a collection of probabilities for a identical node.
    Whenever given r.v.(Random Variable) list, list is regarded Conditional nodes and it return probabilities.

    Example:
      when cpt_type::condition_type cond = {{node1, 2}, {node2, 5}};
           cpt_type cpt(node3, {node1, node2});

      and  std::vector<double> probabilities = cpt[cond];

      then P(node3=0 | node1=2, node2=5): probabilities[0]
           P(node3=1 | node1=2, node2=5): probabilities[1]
                :                                :
           P(node3=n | node1=2, node2=5): probabilities[n]
**/
class cpt {
public:
    using rv_ptr = component::random_variable_ptr;
    using matrix_type = matrix<std::vector<double>>;
    using condition_type = std::unordered_map<rv_ptr, std::size_t>;

    //! (ctor) Initialize CPT by random variable.
    /*! Initialization is performed for only argument random variable.
        The random variable cannot be changed.

        @param[in]  rv: a random variable which represents a target of generated instance.
    **/
    explicit cpt(rv_ptr const& rv)
        : rv_(rv)
    {
        create_matrix({});
    }

    //! (ctor) Initialize CPT by node and parent (conditional) node list
    /*! Initialization is performed for only argument random variab which cannot be changed,
        and generate a condition matrix corresponding to parent nodes in arguments.

        @param[in]  rv: a random variable which represents a target of generated instance.
        @param[in]  parants: a set of random variables corresponding to parent (conditional) nodes of the first argument rv.
    **/
    cpt(rv_ptr const& rv, std::vector<rv_ptr> const& parents)
        : rv_(rv), parents_(parents)
    {
        create_matrix(parents);
    }

    //! (Copy ctor)
    cpt(cpt const& other)
        : rv_(other.rv_), parents_(other.parents_), matrix_(new matrix_type(*other.matrix_))
    {
    }

    //! (Default dtor)
    virtual ~cpt()
    {
        delete matrix_;
    }

    //! (Copy assignment)
    cpt& operator=(cpt const& other)
    {
        *this = cpt{other};
        return *this;
    }

    //! Change parent (conditional) nodes and reset
    /*! Strong guarantee: if an exception is thrown, there are no changes in the CPT.
    */
    void reset(std::vector<rv_ptr> parents)
    {
        create_matrix(parents);
        std::swap(this->parents_, parents_);
    }

    //! Access to probability list by condition r.v. list.
    /*! if condition is invalid, then throw std::runtime_error.

        @return an 1-dimensional probability array corresponding to P(X_Q | X_E1 = x_e1, ..., X_EN = x_en).
    **/
    matrix_type::element_type& at(condition_type const& condition)
    {
        return matrix_->at(condition_to_index(condition));
    }

    //! Access to probability list by condition r.v. list (const version).
    /*! if condition is invalid then throw std::runtime_error.

        @return an 1-dimensional probability array corresponding to P(X_Q | X_E1 = x_e1, ..., X_EN = x_en).
    **/
    matrix_type::element_type const& at(condition_type const& condition) const
    {
        return matrix_->at(condition_to_index(condition));
    }

    //! Access to probability list by condition r.v. list.
    /*! if condition is invalid, then throw std::runtime_error.

        @return an 1-dimensional probability array corresponding to P(X_Q | X_E1 = x_e1, ..., X_EN = x_en).
    **/
    matrix_type::element_type& operator[] (condition_type const& condition)
    {
        return matrix_->operator[](condition_to_index(condition));
    }

    //! Access to probability list by condition r.v. list (const version).
    /*! if condition is invalid then throw std::runtime_error.

        @return an 1-dimensional probability array corresponding to P(X_Q | X_E1 = x_e1, ..., X_EN = x_en).
    **/
    matrix_type::element_type const& operator[] (condition_type const& condition) const
    {
        return matrix_->operator[](condition_to_index(condition));
    }

    //! Get the target random variable of this instance.
    /*! 
        @return a random variable corresponding to this cpt instance.
    **/
    rv_ptr const& rv() const noexcept
    {
        return rv_;
    }
	
	//! Get the parent random variables of the target.
	/*!
		@return a array of random variables which are parents in this cpt instance.
	**/
	std::vector<rv_ptr> const& parents() const noexcept
	{
		return parents_;
	}

	//! Verify whether given condition is sufficient in this instance.
	/*!
		@param[in]   condition: a set of pairs corresponding to the observed (evidence) nodes (environment).
		@return      true if condition is valid for this cpt instance, false if index is invalid.
	**/
	bool is_valid(condition_type const& condition) const noexcept
	{
		for(auto const& parent : parents())
		{
			auto const it = condition.find(parent);
			if(it == condition.end() || it->second >= parent->max_value)
				return false;
		}

		return true;
	}

    //! Convert condition of conditional r.v. into access list of matrix.
    /*! std::out_of_range is thown if the condition does not contain any parent nodes' value.
        std::runtime_error is thown if any parent nodes' value which is contained in the condition exceeds its max value.
		No exceptions are thrown if is_valid(condition) is true.

        @return an array of index in implementation of class cpt.
    **/
    std::vector<std::size_t> condition_to_index(condition_type const& condition) const
    {
        if(parents_.empty())
        {
            return std::vector<std::size_t>{0};
        }
        else
        {
            std::vector<std::size_t> index(parents_.size());
            for(std::size_t i = 0; i < parents_.size(); ++i)
            {
                index[i] = condition.at(parents_[i]); // if condition does not have parents_[i] as key, the exception is thrown.
                if(index[i] >= parents_[i]->max_value)
                    throw std::runtime_error("Condition is out of range");
            }

            return index;
        }
    }

private:
    //! Create matrix<std::vector<double>> table from rv_ and parent_rvs.
    /*! The same matrix's element shows the same condition.
        std::vector<double> show probability that is node correspondence.
        Strong guarantee: if an exception is thrown, there are no changes in the CPT.
    **/
    void create_matrix(std::vector<rv_ptr> const& parents)
    {
        if(rv_->max_value == 0)
            throw std::runtime_error("Target random variables cannot take any value.");

        if(std::any_of(parents.cbegin(), parents.cend(), [](auto const& rv){ return rv->max_value == 0; }))
            throw std::runtime_error("Parent random variables which cannot take any value exist.");

        if(parents.empty())
        {
            auto* new_matrix = new matrix_type(std::vector<double>{1}, std::vector<double>(rv_->max_value));
            std::swap(matrix_, new_matrix);
            delete new_matrix;

        }
        else
        {
            std::vector<std::size_t> size(parents.size());
            std::transform(
                parents.cbegin(), parents.cend(), size.begin(),
                [](rv_ptr const& parent){ return parent->max_value; });

            auto* new_matrix = new matrix_type(size, std::vector<double>(rv_->max_value));
            std::swap(matrix_, new_matrix);
            delete new_matrix;
        }
    }

    rv_ptr const rv_;
    std::vector<rv_ptr> parents_;
    matrix_type* matrix_ = nullptr;
};


//! A manager class of some CPTs corresponding to each nodes.
/*! CPT is unique to each node. One node correspond to one CPT.
**/
class cpt_manager {
public:
    using cpt_type = cpt;
    using node_type = component::node_ptr;

    //! (Default ctor)
    cpt_manager() = default;

    //! (Default dtor)
    virtual ~cpt_manager() = default;

    //! Register CPT to this CPT manager.
    /*! A node corresponding to CPT is given from argument menber function(cpt.node()).

        @return instance of cpt resistered to CPT manager.
    **/
    cpt_type& enroll(node_type const& node, cpt_type const& cpt)
    {
        return enroll(node, cpt_type(cpt));
    }

    //! Register CPT to cpt_list.
    /*! A node corresponding to CPT is given from argument menber function(cpt.node()).
        After call, argument will be moved and invalid.

        @return instance of cpt resistered to CPT manager.
    **/
    cpt_type& enroll(node_type const& node, cpt_type&& cpt)
    {
        auto it = cpt_list_.find(node);
        if(it != cpt_list_.end())
        {
            // override
            it->second = std::move(cpt);
        }
        else
        {
            // add to list, and update iterator variable
            it = cpt_list_.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(node),
                std::forward_as_tuple(std::move(cpt))
                ).first;
        }

        return it->second;
    }

    //! Remove CPT corresponded to argument node from cpt_list.
    /*! If CPT is not found, this function is no effect.
    **/
    void unenroll(node_type const& node)
    {
        auto const it = cpt_list_.find(node);
        if(it != cpt_list_.end())
        {
            cpt_list_.erase(it);
        }
    }

    //! CPT access operator
    /*! If list does not have it, create and resister new CPT.
    **/
    cpt_type& operator[](node_type const& node)
    {
        auto it = cpt_list_.find(node);
        if(it == cpt_list_.end())
        {
            cpt_type new_cpt(node->get());
            it = cpt_list_.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(node),
                std::forward_as_tuple(new_cpt)
                ).first;
        }

        return it->second;
    }

    //! CPT access
    /*! If list does not have it, throw exception.
    **/
    cpt_type& at(node_type const& node)
    {
        return cpt_list_.at(node);
    }

    //! CPT access (const version)
    /*! If list does not have it, throw exception.
    **/
    cpt_type const& at(node_type const& node) const
    {
        return cpt_list_.at(node);
    }

private:
    std::unordered_map<node_type, cpt_type> cpt_list_;
};

} // namespace bn

#endif // #ifndef BAYESIAN_NETWORKS_CPT_HPP
