#ifndef BNI_BP_HPP
#define BNI_BP_HPP

#include <functional>
#include "graph.hpp"

namespace bn {

class bp {
public:
    typedef std::unordered_map<vertex_type, matrix_type> return_type;

    explicit bp(graph_t const& graph);
    virtual ~bp() = default;

    void initialize();

    inline return_type operator()()
    {
        std::unordered_map<vertex_type, matrix_type> const precondition;
        return operator()(precondition);
    }
    return_type operator()(std::unordered_map<vertex_type, matrix_type> const& precondition);

private:
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

    bool is_preconditional_node(vertex_type const& node) const;

    graph_t const graph_;
    std::unordered_map<vertex_type, matrix_type> pi_;
    std::unordered_map<vertex_type, matrix_type> lambda_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> pi_i_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> lambda_k_;
    std::vector<vertex_type> preconditional_node_;

    std::unordered_map<vertex_type, matrix_type> new_pi_;
    std::unordered_map<vertex_type, matrix_type> new_lambda_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> new_pi_i_;
    std::unordered_map<vertex_type, std::unordered_map<vertex_type, matrix_type>> new_lambda_k_;
};

} // namespace bn

#endif // #ifndef BNI_BP_HPP

