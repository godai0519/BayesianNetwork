#ifndef BNI_BP_HPP
#define BNI_BP_HPP

#include <functional>
#include "graph.hpp"

namespace bn {

// 辺のlikelihoodを保持し，簡易に前の状態をコピーすることが可能にするため
class likelihood_list {
public:
    class value_type {
    public:
        value_type(vertex_type const& from, vertex_type const& to, edge_type const& edge);
        virtual ~value_type() = default;

        vertex_type const from() const;
        vertex_type const to() const;
        edge_type const edge() const;

        matrix_type& likelihood();
        matrix_type const& likelihood() const;

    private:
        vertex_type from_;
        vertex_type to_;
        edge_type edge_;
        matrix_type likelihood_;
    };

    matrix_type& add_manage(vertex_type const& from, vertex_type const& to, edge_type const& edge);
    void del_manage(edge_type const& edge);
    void del_manage(vertex_type const& from, vertex_type const& to);

    matrix_type& operator() (edge_type const& edge);
    matrix_type& operator() (vertex_type const& from, vertex_type const& to);
    matrix_type const& operator() (edge_type const& edge) const;
    matrix_type const& operator() (vertex_type const& from, vertex_type const& to) const;

private:
    std::vector<value_type>::iterator find(edge_type const& edge);
    std::vector<value_type>::iterator find(vertex_type const& from, vertex_type const& to);

    std::vector<value_type>::const_iterator find(edge_type const& edge) const;
    std::vector<value_type>::const_iterator find(vertex_type const& from, vertex_type const& to) const;

    std::vector<value_type> data_;
};

class bp {
public:
    bp(graph_t const& graph);
    virtual ~bp() = default;

    matrix_type operator()(
        vertex_type const& target,
        std::vector<std::pair<vertex_type, int>> const& condition
        );

private:
    void calculate_pi(vertex_type const& target);
    void calculate_pi_i(vertex_type const& from, vertex_type const& target);
    void calculate_lambda(vertex_type const& target);
    void calculate_lambda_k(vertex_type const& from, vertex_type const& target);
    
    // 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
    void all_combination_pattern(
        std::vector<vertex_type> const& combination,
        std::function<void(condition_t const&)> const& function
        );

    // 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
    void all_combination_pattern(
        condition_t& combination,
        condition_t::iterator it,
        std::function<void(condition_t const&)> const& func
        );

    // combinationから必要なノードの選択状態だけ取り出して，条件を更新する
    condition_t update_select_condition(condition_t const& whole_condition, condition_t particle_condition);
    condition_t update_select_condition(condition_t const& whole_condition, std::vector<vertex_type> const& particle);

    // 周辺化
    std::unordered_map<condition_t, double> marginalize(
        std::unordered_map<condition_t, double> base,
        vertex_type const& target
        );

    // 条件化(条件を添加する)
    void conditioning(
        std::unordered_map<condition_t, double> const& probabilities,
        vertex_type const& target_node,
        vertex_type const& condition_node
        );

    // 指定ノードから上流に拡散，上流が確定した後に自身を算出することで解を得る
    std::unordered_map<condition_t, double> calculate_likelihood_from_backward(
        graph_t const& graph,
        vertex_type const& node
        );

    // cptを元に全てのエッジのlikelihoodを算出する
    void calculate_likelihood(graph_t const& graph);

    matrix_type propagate_forward(
        graph_t const& graph,
        vertex_type const& target,
        std::vector<std::pair<vertex_type, int>> const& condition
        );
    matrix_type propagate_backward(
        graph_t const& graph,
        vertex_type const& target,
        std::vector<std::pair<vertex_type, int>> const& condition
        );

    graph_t const graph_;
    likelihood_list likelihood_list_;
};

} // namespace bn

#endif // #ifndef BNI_BP_HPP

