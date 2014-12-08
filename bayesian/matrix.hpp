#ifndef BNI_MATRIX_HPP
#define BNI_MATRIX_HPP

#include <vector>
#include <cassert>

namespace bn {

// 制約: 一般的な2次元配列より，アクセスは [y][x]
class matrix_type {
public:
    explicit matrix_type() = default;
    matrix_type(std::size_t const height, std::size_t const width, double const default_value = 0.0)
    {
        resize(height, width, default_value);
    }
    virtual ~matrix_type() = default;

    std::size_t height() const { return height_; }
    std::size_t width() const { return width_; }

    // 指定サイズに伸縮
    void resize(std::size_t const height, std::size_t const width, double const default_value = 0.0)
    {
        // width
        for(auto& line : mat_)
        {
            line.resize(width, default_value);
        }

        // height
        mat_.resize(height, std::vector<double>(width, default_value));

        height_ = height;
        width_  = width;
    }

    // assignに失敗すれば，falseが帰る
    template<class InputIterator>
    bool assign(InputIterator begin, InputIterator const& end)
    {
        if(static_cast<std::size_t>(std::distance(begin, end)) >= width_ * height_)
        {
            for(auto outer_it = mat_.begin(); outer_it != mat_.end(); ++outer_it)
            {
                for(auto inner_it = outer_it->begin(); inner_it != outer_it->end();)
                {
                    *inner_it++ = *begin++;
                }
            }

            return true;
        }
        else
        {
            return false;
        }
    }

    // 内側の配列を返す
    // 設計が雑なの
    std::vector<double>& operator[] (std::size_t const index)
    {
        return mat_[index];
    }
    std::vector<double> const& operator[] (std::size_t const index) const
    {
        return mat_[index];
    }

private:
    std::size_t height_ = 0;
    std::size_t width_ = 0;

    std::vector<std::vector<double>> mat_;
};

} // namespace bn

// 行列同士・行列と整数の積(定義)
bn::matrix_type operator*(bn::matrix_type const& lhs, bn::matrix_type const& rhs);
template<class T> bn::matrix_type operator*(bn::matrix_type const& rhs, T const& lhs);
template<class T> bn::matrix_type operator*(T const& lhs, bn::matrix_type const& rhs);

#include "matrix_impl.hpp"

#endif // #ifndef BNI_MATRIX_HPP
