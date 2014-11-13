#ifndef BNI_MATRIX_HPP
#define BNI_MATRIX_HPP

namespace bn {

template<class T>
class matrix {
public:
    template<class U>
    explicit matrix(U && mat)
      : mat_(std::forward<U>(mat))
    {
    }

    virtual ~matrix() = default;

private:
    std::vector<std::vector<T>> mat_;
};

}

#endif // #ifndef BNI_MATRIX_HPP

