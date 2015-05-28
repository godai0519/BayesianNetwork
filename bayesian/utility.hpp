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

} // namespace bn

#endif // BNI_LEARNING_UTILITY_HPP
