#ifndef BNI_SAMPLING_HPP
#define BNI_SAMPLING_HPP

#include <random>
#include "graph.hpp"

namespace bn {

class sampling {
public:
    typedef std::unordered_map<vertex_type, matrix_type> return_type;
    typedef std::vector<std::unordered_map<vertex_type, int>> pattern_list;

    explicit sampling(graph_t const& graph);
    virtual ~sampling() = default;

    // By-pass
    inline return_type operator()(int const generate_sample_num = 10000)
    {
        std::vector<std::pair<vertex_type, int>> const precondition;
        return operator()(precondition, generate_sample_num);
    }

    // Run: Logic Sampling (a.k.a. Rejection Sampling)
    return_type operator()(
        std::vector<std::pair<vertex_type, int>> const& precondition,
        int const generate_sample_num = 10000
        );

private:
    // 確率に任せたランダム系列のリストを返す
    pattern_list generate_pattern(
        int const num,
        std::vector<std::pair<vertex_type, int>> const& condition = {}
        );

    // 再帰を元にランダムな系列を作成する
    void choice_pattern(
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

    graph_t const graph_;
    probability_generator probability_generator_;
};

} // namespace bn

#endif // #ifndef BNI_SAMPLING_HPP

