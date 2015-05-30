#ifndef BNI_UTILITY_HPP
#define BNI_UTILITY_HPP

#include <string>
#include <random>
#include <vector>
#include <algorithm>
#include <cstdint>

namespace bn {

template<class RandomEngine>
RandomEngine make_engine()
{
    std::random_device rd;

    std::vector<uint32_t> seed_vector(10);
    std::generate(seed_vector.begin(), seed_vector.end(), std::ref(rd));

    std::seed_seq seed(seed_vector.begin(), seed_vector.end());
    return RandomEngine(seed);
}
/*
void all_combination_pattern(
    std::vector<bn::vertex_type> const& combination,
    std::function<void(bn::condition_t const&)> const& function
    )
{
    typedef std::vector<bn::vertex_type>::const_iterator iterator_type;
    std::function<void(iterator_type const, iterator_type const&)> recursive;

    bn::condition_t condition;
    recursive = [&](iterator_type const it, iterator_type const& end)
    {
        if(it == end)
        {
            function(condition);
        }
        else
        {
            for(int i = 0; i < (*it)->selectable_num; ++i)
            {
                condition[*it] = i;
                recursive(it + 1, end);
            }
        }
    };

    recursive(combination.cbegin(), combination.cend());
}
*/

} // namespace bn

#endif // BNI_LEARNING_UTILITY_HPP
