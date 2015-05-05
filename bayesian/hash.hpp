#ifndef BNI_HASH_HPP
#define BNI_HASH_HPP

#include <utility>
#include <bayesian/graph.hpp>

namespace bn {

template<class Key>
class hash;

// For std::unordered_map<std::pair<vertex_type, condition_t>, std::vector<std::size_t>>
template <>
class hash<std::pair<vertex_type, condition_t>> {
public:
    size_t operator()(std::pair<vertex_type, condition_t> const& key) const
    {
        return std::hash<vertex_type>()(key.first) ^ std::hash<condition_t>()(key.second);
    }
};

// For std::unordered_map<std::pair<vertex_type, vertex_type>, double>
template <>
class hash<std::pair<vertex_type, vertex_type>> {
public:
    size_t operator()(std::pair<vertex_type, vertex_type> const& key) const
    {
        return std::hash<vertex_type>()(key.first) ^ std::hash<vertex_type>()(key.second);
    }
};

}

#endif // #ifndef BNI_HASH_HPP
