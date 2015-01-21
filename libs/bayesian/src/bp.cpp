#include <functional>
#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"
#include "bayesian/cpt.hpp" // unordered_mapのネストの解決

namespace bn {

bp::bp(graph_t const& graph)
    : graph_(graph)
{
}

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
                    value *= pi_i[target][xi][0][condition.at(xi)];
                }
                matrix[0][i] += value;
            }
        });

    // 計算された値でpiを更新させる
    pi[target] = matrix;
}

void bp::calculate_pi_i(vertex_type const& from, vertex_type const& target)
{
    auto out_vertexs = graph_.out_vertexs(target);
    out_vertexs.erase(std::find(out_vertexs.cbegin(), out_vertexs.cend(), from));

    // 事前にpiやlambda_kを更新させる
    calculate_pi(target);
    for(auto const& xj : out_vertexs) calculate_lambda_k(xj, target);

    matrix_type matrix = pi[target];
    for(int i = 0; i < target->selectable_num; ++i)
    {
        for(auto const& xj : out_vertexs)
        {
            matrix[0][i] *= lambda_k[xj][target][0][i];
        }
    }

    // 計算された値でpi_iを更新させる
    pi_i[from][target] = matrix;
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
            matrix[0][i] *= lambda_k[xi][target][0][i];
        }
    }

    // 計算された値でlambdaを更新させる
    lambda[target] = matrix;
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
        auto const times = lambda[from][0][i];
        all_combination_pattern(
            in_vertexs,
            [&](condition_t const& cond)
            {
                double value = times * from->cpt[cond].second[i];
                for(auto const& p : cond)
                {
                    if(p.first != target)
                    {
                        value *= pi_i[from][p.first][0][p.second];
                    }
                }
                matrix[0][cond.at(target)] += value;
            });
    }

    // 計算された値でlambdaを更新させる
    lambda_k[from][target] = matrix;
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

// combinationから必要なノードの選択状態だけ取り出して，条件を更新する
condition_t bp::update_select_condition(
    condition_t const& whole_condition,
    condition_t particle_condition
    )
{
    for(auto it = particle_condition.begin(); it != particle_condition.end(); ++it)
    {
        it->second = whole_condition.at(it->first);
    }
    return particle_condition;
}

// 実際こっちであるべき
condition_t bp::update_select_condition(
    condition_t const& whole_condition,
    std::vector<vertex_type> const& particle
    )
{
    condition_t particle_condition;
    for(auto const& v : particle) particle_condition[v] = 0;

    return update_select_condition(whole_condition, particle_condition);

}

// 周辺化
std::unordered_map<condition_t, double> bp::marginalize(
    std::unordered_map<condition_t, double> base,
    vertex_type const& target
    )
{
    std::unordered_map<condition_t, double> result;
    auto condition = base.begin()->first;
    all_combination_pattern(
        condition, condition.begin(),
        [&result, &base, &target](condition_t const& condition)
        {
            auto new_condition = condition;
            new_condition.erase(target);

            if(result.count(new_condition))
            {
                result[new_condition] += base[condition];
            }
            else
            {
                result[new_condition] = base[condition];
            }
        });

    return result;
}


// 条件化(条件を添加する)
void bp::conditioning(
    std::unordered_map<condition_t, double> const& probabilities,
    vertex_type const& target_node,
    vertex_type const& condition_node
    )
{
    auto pair_probabilities = probabilities;
    for(auto it = probabilities.begin()->first.begin(); it != probabilities.begin()->first.end(); ++it)
    {
        if(it->first != target_node && it->first != condition_node)
        {
            pair_probabilities = marginalize(pair_probabilities, it->first);
        }
    }

    for(int i = 0; i < condition_node->selectable_num; ++i)
    {
        double denominator = 0;
        for(int j = 0; j < target_node->selectable_num; ++j)
        {
            auto const target_prob = pair_probabilities[{{condition_node,i},{target_node,j}}];
            denominator += target_prob;
            likelihood_list_(condition_node, target_node)[i][j] = target_prob;
        }

        for(int j = 0; j < target_node->selectable_num; ++j)
        {
            likelihood_list_(condition_node, target_node)[i][j] /= denominator;
        }
    }
}

// 指定ノードから上流に拡散，上流が確定した後に自身を算出することで解を得る
std::unordered_map<condition_t, double> bp::calculate_likelihood_from_backward(
    graph_t const& graph,
    vertex_type const& node
    )
{
    std::vector<vertex_type> direct_parents; // 直接の親
    std::vector<vertex_type> indirect_parents; // 間接の親

    std::vector<std::unordered_map<condition_t, double>> upward_probabilities;

    // 上流ノードの列挙(boost::transformを使うと簡潔)
    for(auto const& upward_edge : graph.in_edges(node))
    {
        auto const upward_node = graph.source(upward_edge);
        direct_parents.push_back(upward_node);

        likelihood_list_.add_manage(upward_node, node, upward_edge);
    }

    // 上流ノードがあるならば，上流を先に確定させる
    // なんでこんなの書いてるんだろうかとイライラしてきた
    for(auto const& upward_node : direct_parents)
    {
        auto const result = calculate_likelihood_from_backward(graph, upward_node);
        if(result.size() == 0) throw std::runtime_error("calculate_likelihood_from_backward: size = 0");
        for(auto const& probability : result.begin()->first)
        {
            // 間接の親かどうか
            auto const& parent_node = probability.first;
            if(std::find(direct_parents.cbegin(), direct_parents.cend(), parent_node) == direct_parents.cend() &&
               std::find(indirect_parents.cbegin(), indirect_parents.cend(), parent_node) == indirect_parents.cend())
            {
                indirect_parents.push_back(parent_node);
            }
        }

        upward_probabilities.push_back(std::move(result));
    }

    // 上流ノード全ての組み合わせを作製して回す
    condition_t combination = {{node, 0}};
    for(auto const& key : direct_parents) combination[key] = 0;
    for(auto const& key : indirect_parents) combination[key] = 0;

    // returnにも使われる同時確率を計算
    std::unordered_map<condition_t, double> target_node_probability;
    all_combination_pattern(
        combination, combination.begin(),
        [this, &node, &target_node_probability, &upward_probabilities](condition_t const& combination)
        {
            // foldl使うと，どうせVC落ちるからやめた
            double probability = node->cpt[update_select_condition(combination, node->cpt.condition_node())]
                                   .second[combination.at(node)];

            for(auto const& upward : upward_probabilities)
            {
               probability *= upward.at(update_select_condition(combination, upward.begin()->first));
            }
            target_node_probability[combination] = probability;
        });

    // 自身と直接の親以外の要素について，周辺化
    auto marginalized_probability = target_node_probability;
    for(auto it = indirect_parents.begin(); it != indirect_parents.end(); ++it)
    {
        marginalized_probability = marginalize(marginalized_probability, *it);
    }

    for(auto it = direct_parents.begin(); it != direct_parents.end(); ++it)
    {
        conditioning(marginalized_probability, node, *it);
    }

    return target_node_probability;
}

// cptを元に全てのエッジのlikelihoodを算出する
void bp::calculate_likelihood(graph_t const& graph)
{
    auto node_list = graph.vertex_list();
    for(auto const& node : node_list)
    {
        auto edges = graph.out_edges(node);
        if(edges.size() == 0)
        {
            // 末端から走査したほうが都合がいいので，末端のノードを全網羅(==全エッジ走査)
            calculate_likelihood_from_backward(graph, node);
        }
    }
}

matrix_type bp::operator()(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    // likelihoodを求める
    calculate_likelihood(graph);

    // 前後の要素に伝播させる
    auto const e_minus = propagate_forward(graph, node, condition);
    auto const e_plus = propagate_backward(graph, node, condition);

    // 掛け算
    auto const elem_num = node->selectable_num;
    double sum = 0.0;
    matrix_type mat(elem_num, 1);
    for(std::size_t i = 0; i < e_minus.height(); ++i)
    {
        double const product = e_minus[i][0] * e_plus[0][i];
        sum += product;
        mat[i][0] = product;
    }

    // 正規化
    for(std::size_t i = 0; i < e_minus.height(); ++i)
    {
        mat[i][0] /= sum;
    }

    return mat;
}

std::pair<bool, int> find_condition(
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    for(auto const& c : condition)
    {
        if(c.first == node)
        {
            return std::make_pair(true, c.second);
        }
    }

    return std::make_pair(false, 0);
}

// 下流要素の確率推論
matrix_type bp::propagate_forward(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    auto const elem_num = node->selectable_num;
    matrix_type mat(elem_num, 1, 1);

    // node ∈ condition
    auto is_condition = find_condition(node, condition);
    if(is_condition.first)
    {
        for(int i = 0; i < elem_num; ++i)
        {
            mat[i][0] = (i == is_condition.second) ? 1 : 0;
        }
        return mat;
    }

    // conditionに含まれないから伝播 (e-要素)
    auto const out_edges = graph.out_edges(node);
    if(!out_edges.empty())
    {
        for(auto const& edge : out_edges)
        {
            mat = mat % (likelihood_list_(edge) * propagate_forward(graph, graph.target(edge), condition));
        }

        return mat;
    }

    // 末端は全ての確率が等しいとする
    for(int i = 0; i < elem_num; ++i)
    {
        mat[i][0] = 1.0;
    }
    return mat;
}

// 上流要素の確率推論
matrix_type bp::propagate_backward(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    auto const elem_num = node->selectable_num;
    matrix_type mat(1, elem_num, 1);

    // node ∈ condition
    auto is_condition = find_condition(node, condition);
    if(is_condition.first)
    {
        for(int i = 0; i < elem_num; ++i)
        {
            mat[0][i] = (i == is_condition.second) ? 1 : 0;
        }
        return mat;
    }

    // conditionに含まれないから伝播 (e+要素)
    auto const in_edges = graph.in_edges(node);
    if(!in_edges.empty())
    {
        for(auto const& edge : in_edges)
        {
            mat = mat % (propagate_backward(graph, graph.source(edge), condition) * likelihood_list_(edge));
        }

        return mat;
    }

    // 最上位ノードは事前確率を割り当てる
    matrix_type evidence(1, node->selectable_num);
    for(int i = 0; i < node->selectable_num; ++i)
    {
		condition_t const cond;
        evidence[0][i] = node->cpt[cond].second[i];
    }
    return evidence;
}

} // namespace bn

