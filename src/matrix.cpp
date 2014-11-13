#include <vector>
#include "matrix.hpp"

namespace bn {

matrix_type operator*(matrix_type const& lhs, matrix_type const& rhs)
{
    assert(lhs.shape()[1] == rhs.shape()[0]);

    matrix_type result(boost::extents[lhs.shape()[0]][rhs.shape()[1]]); // 0.0 fill (http://melpon.org/wandbox/permlink/jXd8ETmx8y8ZwK4A)
    for(matrix_type::size_type i = 0; i < result.shape()[0]; ++i)
    {
        for(matrix_type::size_type j = 0; j < result.shape()[1]; ++j)
        {
            for(matrix_type::size_type k = 0; k < lhs.shape()[1]; ++k)
            {
                result[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }

    return result;
}

} // namespace bn
