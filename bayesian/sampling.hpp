#ifndef BNI_SAMPLING_HPP
#define BNI_SAMPLING_HPP

#include <random>
#include "graph.hpp"

namespace bn {

class sampling {
public:
    typedef std::unordered_map<vertex_type, matrix_type> return_type;
    typedef std::vector<std::unordered_map<vertex_type, int>> pattern_list;

    sampling() = default;
    virtual ~sampling() = default;

    return_type operator()(
        graph_t const& graph,
        std::vector<std::pair<vertex_type, int>> const& condition
        );

private:
    // 確率に任せたランダム系列のリストを返す
    pattern_list generate_pattern(
        graph_t const& graph,
        int const num,
        std::vector<std::pair<vertex_type, int>> const& condition = {}
        );

    // 再帰を元にランダムな系列を作成する
    void choice_pattern(
        graph_t const& graph,
        vertex_type const& current,
        std::vector<vertex_type>& remain_node,
        std::unordered_map<vertex_type, int>& pattern
        );

    // 乱数生成器
    class probability_generator {
    public:
        probability_generator();
        double operator() ();
    private:
        std::unique_ptr<std::mt19937> engine_;
        std::uniform_real_distribution<double> distribution_;
    };

    probability_generator probability_generator_;
};

} // namespace bn

#endif // #ifndef BNI_SAMPLING_HPP

