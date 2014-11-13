#include "matrix.hpp"

int main()
{
    bn::matrix_type mat1(boost::extents[2][3]);
    bn::matrix_type mat2(boost::extents[3][2]);

    mat1[0][0] = 1;
    mat1[0][1] = 3;
    mat1[0][2] = 5;
    mat1[1][0] = 2;
    mat1[1][1] = 4;
    mat1[1][2] = 6;

    mat2[0][0] = 1;
    mat2[0][1] = 2;
    mat2[1][0] = 3;
    mat2[1][1] = 4;
    mat2[2][0] = 5;
    mat2[2][1] = 6;

    auto res = 10 * mat2;
    for(bn::matrix_type::size_type i = 0; i < res.shape()[0]; ++i)
    {
        for(bn::matrix_type::size_type j = 0; j < res.shape()[1]; ++j)
        {
            std::cout << res[i][j] << " ";
        }
        std::cout << "\n";
    }

}

