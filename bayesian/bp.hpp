#ifndef BNI_BP_HPP
#define BNI_BP_HPP

#include <functional>
#include "graph.hpp"

namespace bn {

class bp {
public:
    typedef std::unordered_map<vertex_type, matrix_type> return_type;

    bp(graph_t const& graph);
    virtual ~bp() = default;

    return_type operator()(
/*
        vertex_type const& target,
        std::vector<std::pair<vertex_type, int>> const& condition
*/
        );

private:
    void initialize();
    matrix_type& normalize(matrix_type& target);

    void calculate_pi(vertex_type const& target);
    void calculate_pi_i(vertex_type const& from, vertex_type const& target);
    void calculate_lambda(vertex_type const& target);
    void calculate_lambda_k(vertex_type const& from, vertex_type const& target);
    
    // 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
    void all_combination_pattern(
        std::vector<vertex_type> const& combination,
        std::function<void(condition_t const&)> const& function
        );

    graph_t const graph_;
    std::unordered_map<vertex_type, matrix_type> pi_;
    std::unordered_map<vertex_type, matrix_type> lambda_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> pi_i_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> lambda_k_;
};

} // namespace bn

#endif // #ifndef BNI_BP_HPP

