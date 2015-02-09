#ifndef BNI_BP_HPP
#define BNI_BP_HPP

#include <functional>
#include "graph.hpp"

namespace bn {

class belief_propagation {
public:
    typedef std::unordered_map<vertex_type, matrix_type> return_type;

    explicit belief_propagation(graph_t const& graph);
    virtual ~belief_propagation() = default;

    // By-pass
    inline return_type operator()(double const epsilon = 0.001)
    {
        std::unordered_map<vertex_type, matrix_type> const precondition;
        return operator()(precondition, epsilon);
    }

    // Run: Loopy Belief Propagation
    return_type operator()(std::unordered_map<vertex_type, matrix_type> const& precondition, double const epsilon = 0.001);

private:
    void initialize();
    void calculate_pi(vertex_type const& target);
    void calculate_pi_i(vertex_type const& from, vertex_type const& target);
    void calculate_lambda(vertex_type const& target);
    void calculate_lambda_k(vertex_type const& from, vertex_type const& target);

    // 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
    void all_combination_pattern(
        std::vector<vertex_type> const& combination,
        std::function<void(condition_t const&)> const& function
        ) const;

    // matrixに存在するすべての数の和が1になるように正規化する
    matrix_type normalize(matrix_type target) const;

    // 引数のノードが事前条件ノードであるかどうか確認する
    bool is_preconditional_node(vertex_type const& node) const;

    graph_t const graph_;
    std::vector<vertex_type> preconditional_node_;

    // 現在のメッセージ状況
    std::unordered_map<vertex_type, matrix_type> pi_;
    std::unordered_map<vertex_type, matrix_type> lambda_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> pi_i_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> lambda_k_;

    // 未来のメッセージ状況
    std::unordered_map<vertex_type, matrix_type> new_pi_;
    std::unordered_map<vertex_type, matrix_type> new_lambda_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> new_pi_i_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> new_lambda_k_;
};

} // namespace bn

#endif // #ifndef BNI_BP_HPP
