/**
* @file matrix.hpp
* @brief Provide a class representing matrix.
* @author godai_0519
* @date 03/05/2017
*/

#ifndef BAYESIAN_NETWORKS_MATRIX_HPP
#define BAYESIAN_NETWORKS_MATRIX_HPP

#include <cstdint>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <functional>

namespace bn {

//! Provide a class representing dinamic N dimensions (N is equal to or greater than 1).
/*! Dimensions can be determined at run time.
    Element access can be performed by std::vector<std::size_t> instead of operetors [][]...[].

    @tparam Elem: Type of elements hold by matrix<Elem>.
**/
template<class Elem>
class matrix {
public:
    typedef Elem element_type;
    typedef std::vector<std::size_t> sizes_type;
    typedef std::vector<std::size_t> steps_type;
    typedef std::vector<element_type> data_type;

    //! (Default ctor) Initialize matrix as empty.
    explicit matrix()
        : sizes_(), capacities_(), steps_(), data_()
    {
    }

    //! (ctor) Initialize matrix according to array of size-hint.
    /*! @tparam T: Type of list containing size information.
        @param[in]   sizes: a list containing sizes of each dimension in matrix,
                            such that, sizes.size() == dimention.
    **/
    template<class T>
    explicit matrix(T const& sizes)
        : sizes_(std::cbegin(sizes), std::cend(sizes))
        , capacities_(std::cbegin(sizes), std::cend(sizes))
        , steps_(calc_steps(capacities_))
        , data_(calc_elem_size(sizes))
    {
        if(sizes_.empty())
            throw std::runtime_error("bn::matrix do not allow zero dimension.")
    }

    //! (ctor) Initialize matrix by filling with default_value according to array of size-hint.
    /*! The same function of explicit matrix(T const& size) exist;
        however, every data is filled by default_value in this ctor.

        @tparam T: Type of list containing size information.
        @param[in]   sizes: a list containing sizes of each dimension in matrix,
                            such that, sizes.size() == dimention.
        @param[in]   default_value: a value used for filling.
    **/
    template<class T>
    matrix(T const& sizes, Elem const& default_value)
        : sizes_(std::cbegin(sizes), std::cend(sizes))
        , capacities_(std::cbegin(sizes), std::cend(sizes))
        , steps_(calc_steps(capacities_))
        , data_(calc_elem_size(sizes), default_value)
    {
        if (sizes_.empty())
            throw std::runtime_error("bn::matrix do not allow zero dimension.")
    }

    //! (ctor) Initialize matrix by filling with iterator according to array of size-hint.
    /*! The same function of explicit matrix(T const& size) exist;
        however, in this ctor, every data is filled by values of range [first, last).

        @tparam T: Type of list containing size information.
        @tparam InputIterator: Type of iterator containing data to fill.
        @param[in]   sizes: a list containing sizes of each dimension in matrix,
                            such that, sizes.size() == dimention.
        @param[in]   first: iterator of begining of data to filling.
        @param[in]   last: iterator of ending of data to filling.
    **/
    template<class T, class InputIterator>
    matrix(T const& sizes, InputIterator first, InputIterator last)
    {
        sizes_.assign(std::cbegin(sizes), std::cend(sizes));
        capacities_.assign(std::cbegin(sizes), std::cend(sizes));
        steps_ = calc_steps(capacities_);
        if (sizes_.empty())
            throw std::runtime_error("bn::matrix do not allow zero dimension.")

        auto const data_size = calc_elem_size(sizes);
        auto const iterator_size = std::distance(first, last);
        if(data_size != iterator_size)
            throw std::runtime_error("matrix::matrix: iterator data size is invalid");

        data_.assign(first, last);
    }

    //! (Copy ctor)
    matrix(matrix<Elem> const&) = default;

    //! (Move ctor)
    matrix(matrix<Elem>&&) = default;

    //! (Default dtor)
    virtual ~matrix() = default;

    //! Assign matrix by filling with default_value according to array of size-hint.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the matrix.

        @tparam T: Type of list containing size information.
        @param[in]   sizes: a list containing sizes of each dimension in matrix,
                            such that, sizes.size() == dimention.
        @param[in]   default_value: a value used for filling.
    **/
    template<class T>
    void assign(T const& sizes, Elem const& default_value)
    {
        *this = matrix<Elem>{sizes, default_value };
    }

    //! Assign matrix by filling with iterator according to array of size-hint.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the matrix.

        @tparam T: Type of list containing size information.
        @tparam InputIterator: Type of iterator containing data to fill.
        @param[in]   sizes: a list containing sizes of each dimension in matrix,
                            such that, sizes.size() == dimention.
        @param[in]   first: iterator of begining of data to filling.
        @param[in]   last: iterator of ending of data to filling.
    **/
    template<class T, class InputIterator>
    void assign(T const& sizes, InputIterator first, InputIterator last)
    {
        *this = matrix<Elem>{sizes, first, last};
    }

    //! (Copy operator=) Initialize nodes and arcs by copying a parameter.
    /*! Strong guarantee: if an exception is thrown, there are no changes in the matrix.
    **/
    matrix<Elem>& operator=(matrix<Elem> const& rhs)
    {
        *this = matrix<Elem>{rhs};
        return *this;
    }

    //! (Move operator=) Initialize nodes and arcs by copying a parameter.
    /*! The parameter will be destroyed.
        Strong guarantee: if an exception is thrown, there are no changes in the matrix.
    **/
    matrix<Elem>& operator=(matrix<Elem>&&) = default;

    //! Clear martix's all data.
    /*! This function invelidate member variables' any past-the-end iterators, pointers, and references
        based on the specification of std::vector.
        No-throw guarantee: never throws exceptions.
    **/
    void clear() noexcept
    {
        sizes_.clear();
        capacities_.clear();
        steps_.clear();
        data_.clear();
    }

    //! Resize this container by argument sizes
    /*! If a element of sizes which is larger than capacity corresponding to the same index is exist,
        matrix is resized by allocating new memory space; otherwise, size_ is updated while keeping capacities_.
        Strong guarantee: if an exception is thrown, there are no changes in the matrix.

        @param[in]   sizes: a list containing sizes of each dimension in matrix.
    **/
    template<class T>
    void resize(T const& sizes)
    {
        // Zero dimension is not allowed
        if(std::distance(std::cbegin(sizes), std::cend(sizes)) == 0)
            throw std::runtime_error("matrix::size: Zero dimension is not allowed.");

        bool is_necessary_realloc = (sizes.size() > capacities_.size()); // not enough dimension of capacity to become large.
        for(std::size_t i = 0; !is_necessary_realloc && i < sizes.size(); ++i)
        {
            // not enough size of each dimension to become large
            if(sizes[i] > capacities_[i])
                is_necessary_realloc = true;
        }

        if(is_necessary_realloc)
        {
            // One or more Dimension
            sizes_type new_sizes(std::cbegin(sizes), std::cend(sizes));
            sizes_type new_capacities(new_sizes);
            steps_type new_steps = calc_steps(new_capacities);

            // Copy elements of old data into new data.
            data_type new_data(calc_elem_size(sizes));
            for_all([&new_data, &new_steps](sizes_type const& index)
            {
                std::size_t const old_index = continuous_index_from(std::cbegin(index), std::cend(index), std::cbegin(steps_));
                std::size_t const new_index = continuous_index_from(std::cbegin(index), std::cend(index), std::cbegin(new_steps));
                new_data[new_index] = data_[old_index];
            });

            sizes_ = std::move(new_sizes);
            capacities_ = std::move(new_capacities);
            steps = std::move(new_steps);
            data_ = std::move(new_data);
        }
        else
        {
            // Reallocatable is not required.
            // Keep the values of capacity, and Reduce the values of size.
            for(std::size_t i = 0; i < sizes.size(); ++i)
                sizes_[i] = sizes[i];

            for(std::size_t i = sizes.size(); i < sizes_.size(); ++i)
                sizes_[i] = 0;
        }
    }

    //! Resize this container by argument sizes when filling extended elements by value
    /*! If a element of sizes which is larger than capacity corresponding to the same index is exist,
        matrix is resized by allocating new memory space; otherwise, size_ is updated while keeping capacities_.
        Strong guarantee: if an exception is thrown, there are no changes in the matrix.

        @param[in]   sizes: a list containing sizes of each dimension in matrix.
        @param[in]   value: a value which is used for filling extended elements.
    **/
    template<class T>
    void resize(T const& sizes, Elem const& value)
    {
        // Zero dimension is not allowed
        if (std::distance(std::cbegin(sizes), std::cend(sizes)) == 0)
            throw std::runtime_error("matrix::size: Zero dimension is not allowed.");

        bool is_necessary_realloc = (sizes.size() > capacities_.size()); // not enough dimension of capacity to become large.
        for (std::size_t i = 0; !is_necessary_realloc && i < sizes.size(); ++i)
        {
            // not enough size of each dimension to become large
            if (sizes[i] > capacities_[i])
                is_necessary_realloc = true;
        }

        if (is_necessary_realloc)
        {
            // One or more Dimension
            sizes_type new_sizes(std::cbegin(sizes), std::cend(sizes));
            sizes_type new_capacities(new_sizes);
            steps_type new_steps = calc_steps(new_capacities);

            // Copy elements of old data into new data.
            data_type new_data(calc_elem_size(sizes), value);
            for_all([&new_data, &new_steps](sizes_type const& index)
            {
                std::size_t const old_index = continuous_index_from(std::cbegin(index), std::cend(index), std::cbegin(steps_));
                std::size_t const new_index = continuous_index_from(std::cbegin(index), std::cend(index), std::cbegin(new_steps));
                new_data[new_index] = data_[old_index];
            });

            sizes_ = std::move(new_sizes);
            capacities_ = std::move(new_capacities);
            steps = std::move(new_steps);
            data_ = std::move(new_data);
        }
        else
        {
            // Reallocatable is not required.
            // Keep the values of capacity, and Reduce the values of size.
            for (std::size_t i = 0; i < sizes.size(); ++i)
                sizes_[i] = sizes[i];

            for (std::size_t i = sizes.size(); i < sizes_.size(); ++i)
                sizes_[i] = 0;
        }
    }

    //! Request the removal of unused capacity.
    /*! This function reduce capacities to size.
        If any variable of Elem can be swapped with another by no-throw guarantee, the function never throws exceptions (no-throw guarantee).
        Otherwise, if an exception is thrown, this instance is in a valid state (Basic guarantee).
    **/
    void shrink_to_fit()
    {
        sizes_type new_capacities(sizes_);
        steps_type new_steps(calc_steps(sizes_));

        for(std::size_t i = 0; index < data_.size(); ++index)
        {
            sizes_type index(sizes_.size());
            index_from(i, index.begin(), index.end());

            auto const new_i = continuous_index_from(index, new_capacities.cbegin(), new_capacities.cend());
            std::swap(data_[i], data_[new_i]);
        }

        capacities_ = std::move(new_capacities);
        steps_ = std::move(new_steps);
    }

    //! Access to element.
    /*! If index is invalid, then throw std::runtime_error.
        @param[in]   index: std::vector, std::array, or more (overloaded std::cbegin, std::cend).
    **/
    template<class T>
    Elem& at(T const& index)
    {
        if(boundary_check(std::cbegin(index), std::cend(index)))
            throw std::runtime_error("bn::matrix: out of bounds"); // TODO: User Difinition Exception

        return operator[](index);
    }

    //! Access to element (const version).
    /*! If index is invalid, then throw std::runtime_error.

        @param[in]   index: std::vector, std::array, or more (overloaded std::cbegin, std::cend).
    **/
    template<class T>
    Elem const& at(T const& index) const
    {
        if(boundary_chack(std::cbegin(index), std::cend(index)))
            throw std::runtime_error("bn::matrix: out of bounds");

        return operator[](index);
    }

    //! Access to element.
    /*! If index is invalid, this function encounter undefined behavior.

        @param[in]   index: std::vector, std::array, or more (overloaded std::cbegin, std::cend).
    **/
    template<class T>
    Elem& operator[] (T const& index)
    {
        return data_[continuous_index_from(std::cbegin(index), std::cend(index))];
    }

    //! Access to element (const version).
    /*! If index is invalid, this function encounter undefined behavior.

        @param[in]   index: std::vector, std::array, or more (overloaded std::cbegin, std::cend).
    **/
    template<class T>
    Elem const& operator[] (T const& index) const
    {
        return data_[continuous_index_from(std::cbegin(index), std::cend(index))];
    }

    //! Get the number of matrix dimentions.
    std::size_t dims() const noexcept
    {
        return sizes_.size();
    }

    //! Get size of each dimentions.
    sizes_type const& sizes() const noexcept
    {
        return sizes_;
    }

    //! Get capacity of each dimentions.
    sizes_type const& capacities() const noexcept
    {
        return capacities_;
    }

    //! Get steps of each dimentions.
    steps_type const& steps() const noexcept
    {
        return steps_;
    }

    //! Access to const raw data.
    data_type const& data() noexcept
    {
        return data_;
    }

    //! Access to const raw data.
    data_type const& data() const noexcept
    {
        return data_;
    }

    //! Calculate matrix inner product (dot).
    /*! This is left hand, other is right hand.

        @param[in]   other: other or this instance of bn::matrix.
    **/
    template<class T>
    double dot(matrix<T> const& other) const
    {
        if(this->dims() != 1 || other.dims() != 1)
            throw std::runtime_error("matrix::dot requires 1 dimension");

        if(this->sizes()[0] != other.sizes()[0])
            throw std::runtime_error("matrix::dot requires equal length");

        return std::inner_product(this->data_.cbegin(), this->data_.cend(), other.data_.cbegin(), 0.0);
    }

private:
    //! Check validity of index (iterator).
    /*! If within the range of validity then false is returned; otherwise, true is returned.

        @param[in]   index_begin: begin of iterator of index container.
        @param[in]   index_end: nd of iterator of index container.
        @return      false if index is valid, true if index is invalid.
    **/
    template<class Iterator>
    bool boundary_check(Iterator index_begin, Iterator const& index_end) const
    {
        auto sizes_it = sizes_.begin();
        while(index_begin != index_end)
        {
            if(*index_begin++ >= *sizes_it++)
                return true;
        }

        return false;
    }

    // calculate linear index by steps_ and index (iterator)
    template<class Iterator>
    std::size_t continuous_index_from(Iterator index_begin, Iterator const& index_end) const
    {
        return std::inner_product(index_begin, index_end, steps_.begin(), (std::size_t)0);
    }

    // calculate linear index by steps_ and index (iterator)
    template<class InputIterator1, class InputIterator2>
    std::size_t continuous_index_from(InputIterator1 index_begin, InputIterator1 index_end, InputIterator2 step_begin) const
    {
        return std::inner_product(index_begin, index_end, step_begin, (std::size_t)0);
    }

    template<class OutputIterator>
    bool index_from(std::size_t index, OutputIterator index_begin, OutputIterator const& index_end) const
    {
        for(std::size_t i = 0; i < steps_.size(); ++i)
        {
            // errored by not enough container size corresponding iterator
            if(index_begin == index_end) return false;

            *index_begin++ = index / steps_[i];
            index          = index % steps_[i];
        }
        return true;
    }

    // calculate steps from size array
    template<class T>
    steps_type calc_steps(T const& sizes) const
    {
        auto sizes_reverse_first = std::crbegin(sizes);
        auto const sizes_reverse_last = std::crend(sizes);
        auto const sizes_length = std::distance(sizes_reverse_first, sizes_reverse_last);

        if(sizes_length == 0)
        {
            return steps_type(1, 0);
        }
        else
        {
            steps_type steps(sizes_length);
            auto steps_it = steps.rbegin();

            std::size_t data_size = 1;
            while(sizes_reverse_first != sizes_reverse_last)
            {
                *steps_it = data_size;
                data_size *= *sizes_reverse_first;

                ++sizes_reverse_first;
                ++steps_it;
            }

            return steps;
        }
    }

    // calculate steps from size array
    template<class T>
    std::size_t calc_elem_size(T const& sizes) const
    {
        return std::accumulate(std::cbegin(sizes), std::cend(sizes), 1, std::multiplies<std::size_t>());
    }

    void for_all(std::function<void(sizes_type)> const& func) const
    {
        sizes_type index(steps_.size());
        for(std::size_t i = 0; i < data_.size(); ++i)
        {
            index_from(i, index.begin(), index.end());
            for(std::size_t j = 0; j < steps_.size(); ++j)
                if(index[j] >= sizes_[j])
                    continue;

            func(index);
        }
    }


    sizes_type sizes_;
    sizes_type capacities_;
    steps_type steps_;
    data_type data_;
};


//! Calculate matrix inner product (dot) (non-member function).
/*! @param[in]   lhs: This is Left hand of dot.
    @param[in]   rhs: This is Right hand of dot.
    @return      result of doting vector.
 **/
template<class T, class U>
double dot(matrix<T> const& lhs, matrix<U> const& rhs)
{
    return lhs.dot(rhs);
}


//! Calculate matrix-multiplication.
/*! Matrix is multiplied by matrix.

    @param[in]   lhs: This is Left hand matrix of multiplication.
    @param[in]   rhs: This is Right hand matrix of multiplication.
    @return      result of multiplication of matrix.
**/
template<class T, class U>
matrix<decltype((*((T*)0)) * (*((U*)0)))> operator* (matrix<T> const& lhs, matrix<U> const& rhs)
{
    if(lhs.dims() != 2 || rhs.dims() != 2)
        throw std::runtime_error("matrix::operator* requires 2 dimension");

    if(lhs.sizes()[1] != rhs.sizes()[0])
        throw std::runtime_error("matrix::operator* requires lhs.sizes()[1] == rhs.sizes()[0]");

    // create a result container
    std::array<std::size_t, 2> result_size{{lhs.sizes()[0], rhs.sizes()[1]}};
    matrix<decltype((*((T*)0)) + (*((U*)0)))> result(result_size);

    // calc
    std::array<std::size_t, 2> index;
    for(index[0] = 0; index[0] < result_size[0]; ++index[0])
    {
        for(index[1] = 0; index[1] < result_size[1]; ++index[1])
        {
            decltype((*((T*)0)) + (*((U*)0))) value = 0;
            for(std::size_t i = 0; i < lhs.sizes()[1]; ++i)
            {
                value += lhs[std::array<std::size_t, 2>{{index[0], i}}] *
                         rhs[std::array<std::size_t, 2>{{i, index[1]}}];
            }
            result[index] = value;
        }
    }

    return result;
}


//! Calculate matrix-multiplication.
/*! Matrix is multiplied by scala.

    @param[in]   lhs: This is Left hand matrix of multiplication.
    @param[in]   rhs: This is Right hand scala of multiplication.
    @return      result of multiplication of matrix and scala.
**/
template<class T, class U>
matrix<decltype((*((T*)0)) * (*((U*)0)))> operator* (matrix<T> const& lhs, U const rhs)
{
    // calc
    auto const& data = lhs.data();
    std::vector<decltype((*((T*)0)) * (*((U*)0)))> result_data(data.size());
    std::transform(
        data.begin(), data.end(), result_data.begin(),
        [rhs](T const value){ return value * rhs; }
        );

    // create a result container
    matrix<decltype((*((T*)0)) * (*((U*)0)))> result(lhs.sizes(), result_data.begin(), result_data.end());
    return result;
}


//! Calculate matrix-multiplication.
/*! Scala is multiplied by matrix.

    @param[in]   lhs: This is Left hand scala of multiplication.
    @param[in]   rhs: This is Right hand matrix of multiplication.
    @return      result of multiplication of matrix and scala.
**/
template<class T, class U>
matrix<decltype((*((T*)0)) * (*((U*)0)))> operator* (T const lhs, matrix<U> const& rhs)
{
    return rhs * lhs;
}

} // namespace bn

#endif // #ifndef BNI_MATRIX_HPP
