#include <functional>
#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"
#include "bayesian/cpt.hpp" // unordered_mapのネストの解決

namespace bn {

bp::bp(graph_t const& graph)
    : graph_(graph)
{
}

/*
matrix_type bp::operator()(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
}
*/

// TODO: 条件の指定されたノードであった場合の処理
void bp::calculate_pi(vertex_type const& target)
{
    auto const in_vertexs = graph_.in_vertexs(target);
    
    // 事前にpi_iを更新させる
    for(auto const& xi : in_vertexs) calculate_pi_i(target, xi);

    matrix_type matrix(1, target->selectable_num, 0.0);
    all_combination_pattern(
        in_vertexs,
        [&](condition_t const& condition)
        {
            auto const& cpt_data = target->cpt[condition].second;
            for(int i = 0; i < target->selectable_num; ++i)
            {
                double value = cpt_data[i];
                for(auto const& xi : in_vertexs)
                {
                    value *= pi_i_[target][xi][0][condition.at(xi)];
                }
                matrix[0][i] += value;
            }
        });

    // 計算された値でpiを更新させる
    pi_[target] = matrix;
}

void bp::calculate_pi_i(vertex_type const& from, vertex_type const& target)
{
    auto out_vertexs = graph_.out_vertexs(target);
    out_vertexs.erase(std::find(out_vertexs.cbegin(), out_vertexs.cend(), from));

    // 事前にpiやlambda_kを更新させる
    calculate_pi(target);
    for(auto const& xj : out_vertexs) calculate_lambda_k(xj, target);

    matrix_type matrix = pi_[target];
    for(int i = 0; i < target->selectable_num; ++i)
    {
        for(auto const& xj : out_vertexs)
        {
            matrix[0][i] *= lambda_k_[xj][target][0][i];
        }
    }

    // 計算された値でpi_iを更新させる
    pi_i_[from][target] = matrix;
}

void bp::calculate_lambda(vertex_type const& target)
{
    auto const out_vertexs = graph_.out_vertexs(target);

    // 事前にlambda_kを更新させる
    for(auto const& xi : out_vertexs) calculate_lambda_k(xi, target);

    matrix_type matrix(1, target->selectable_num, 1.0);
    for(int i = 0; i < target->selectable_num; ++i)
    {
        for(auto const& xi : out_vertexs)
        {
            matrix[0][i] *= lambda_k_[xi][target][0][i];
        }
    }

    // 計算された値でlambdaを更新させる
    lambda_[target] = matrix;
}

void bp::calculate_lambda_k(vertex_type const& from, vertex_type const& target)
{
    auto const in_vertexs = graph_.in_vertexs(from);
    
    // 事前にlambdaを更新させる
    calculate_lambda(from);
    for(auto const& xl : in_vertexs) calculate_pi_i(from, xl);

    matrix_type matrix(1, target->selectable_num, 0.0);
    for(int i = 0; i < from->selectable_num; ++i)
    {
        auto const times = lambda_[from][0][i];
        all_combination_pattern(
            in_vertexs,
            [&](condition_t const& cond)
            {
                double value = times * from->cpt[cond].second[i];
                for(auto const& p : cond)
                {
                    if(p.first != target)
                    {
                        value *= pi_i_[from][p.first][0][p.second];
                    }
                }
                matrix[0][cond.at(target)] += value;
            });
    }

    // 計算された値でlambdaを更新させる
    lambda_k_[from][target] = matrix;
}

// 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
void bp::all_combination_pattern(
    std::vector<vertex_type> const& combination,
    std::function<void(condition_t const&)> const& function
    )
{
    typedef std::vector<vertex_type>::const_iterator iterator_type;
    std::function<void(iterator_type const, iterator_type const&)> recursive;

    condition_t condition;
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

} // namespace bn

