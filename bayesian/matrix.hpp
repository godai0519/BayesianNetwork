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

    // 要素同士の積 result(i,j) = lhs(i,j) * rhs(i,j)
    // 外積(cross)でも内積(dot)でもない
    matrix_type& operator%=(matrix_type const& rhs)
    {
        assert(this->width() == rhs.width() && this->height() == rhs.height());

        for(std::size_t i = 0; i < this->height(); ++i)
        {
            for(std::size_t j = 0; j < this->width(); ++j)
            {
                mat_[i][j] *= rhs.mat_[i][j];
            }
        }

        return *this;
    }

    matrix_type operator%(matrix_type const& rhs) const
    {
        matrix_type tmp(*this);
        tmp %= rhs;
        return tmp;
    }

    // 内積(dot)
    bn::matrix_type operator*=(bn::matrix_type const& rhs)
    {
        assert(this->width() == rhs.height()); // requirement for calculating the product

        std::vector<std::vector<double>> result(this->height(), std::vector<double>(rhs.width(), 0.0));
        for(std::size_t i = 0; i < this->height(); ++i)
        {
            for(std::size_t j = 0; j < rhs.width(); ++j)
            {
                for(std::size_t k = 0; k < rhs.height(); ++k)
                {
                    result[i][j] += this->mat_[i][k] * rhs.mat_[k][j];
                }
            }
        }

        width_ = rhs.width();
        mat_ = std::move(result);

        return *this;
    }

    bn::matrix_type operator*(bn::matrix_type const& rhs) const
    {
        matrix_type tmp(*this);
        tmp *= rhs;
        return tmp;
    }

private:
    std::size_t height_ = 0;
    std::size_t width_ = 0;

    std::vector<std::vector<double>> mat_;
};

} // namespace bn

namespace {
    using bn::matrix_type;
}

// 行列の整数倍
template<class T>
matrix_type operator*(matrix_type const& rhs, T const& lhs)
{
    // after copy, process all elem
    auto result = rhs;
    for(std::size_t i = 0; i < result.height(); ++i)
    {
        for(std::size_t j = 0; j < result.width(); ++j)
        {
            result[i][j] *= lhs;
        }
    }
    return result;
}

// delegate "integer * matrix" to "matrix * integer"
template<class T>
matrix_type operator*(T const& lhs, matrix_type const& rhs)
{
    return rhs * lhs;
}

#endif // #ifndef BNI_MATRIX_HPP
