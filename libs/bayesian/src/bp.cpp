#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"
#include "cpt.hpp" // unordered_mapのネストの解決

namespace bn {

// 辺のlikelihoodを保持し，簡易に前の状態をコピーすることが可能にするため
class likelihood_list {
public:
    class value_type {
    public:
        value_type(vertex_type const& from, vertex_type const& to, edge_type const& edge)
            : from_(from), to_(to), edge_(edge), likelihood_(from->selectable_num, to->selectable_num)
        {
        }
        virtual ~value_type() = default;

        vertex_type const from() const { return from_; }
        vertex_type const to() const { return to_; }
        edge_type const edge() const { return edge_; }

        matrix_type& likelihood()
        {
            return likelihood_;
        }
        matrix_type const& likelihood() const
        {
            return likelihood_;
        }

    private:
        vertex_type from_;
        vertex_type to_;
        edge_type edge_;
        matrix_type likelihood_;
    };

    matrix_type& add_manage(vertex_type const& from, vertex_type const& to, edge_type const& edge)
    {
        data_.emplace_back(from, to, edge);
        return data_.back().likelihood();
    }

    void del_manage(edge_type const& edge)
    {
        auto it = find(edge);
        if(it != data_.end()) data_.erase(it);
    }

    void del_manage(vertex_type const& from, vertex_type const& to)
    {
        auto it = find(from, to);
        if(it != data_.end()) data_.erase(it);
    }

    matrix_type& operator() (edge_type const& edge)
    {
        auto it = find(edge);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (edge_type)");
        return it->likelihood();
    }

    matrix_type const& operator() (edge_type const& edge) const
    {
        auto it = find(edge);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (edge_type) const");
        return it->likelihood();
    }

    matrix_type& operator() (vertex_type const& from, vertex_type const& to)
    {
        auto it = find(from, to);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (vertex_type,vertex_type)");
        return it->likelihood();
    }

    matrix_type const& operator() (vertex_type const& from, vertex_type const& to) const
    {
        auto it = find(from, to);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (vertex_type,vertex_type) const");
        return it->likelihood();
    }

private:
    std::vector<value_type>::iterator find(edge_type const& edge)
    {
        auto it = data_.begin();
        while(it != data_.end())
        {
            if(it->edge() == edge) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type>::iterator find(vertex_type const& from, vertex_type const& to)
    {
        auto it = data_.begin();
        while(it != data_.end())
        {
            if(std::tie(it->from(), it->to()) == std::tie(from, to)) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type>::const_iterator find(edge_type const& edge) const
    {
        auto it = data_.cbegin();
        while(it != data_.cend())
        {
            if(it->edge() == edge) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type>::const_iterator find(vertex_type const& from, vertex_type const& to) const
    {
        auto it = data_.cbegin();
        while(it != data_.cend())
        {
            if(std::tie(it->from(), it->to()) == std::tie(from, to)) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type> data_;
};

// all_combination
void all_combination_pattern(
    std::unordered_map<vertex_type, int>& combination,
    std::unordered_map<vertex_type, int>::iterator it,
    std::function<void(std::unordered_map<vertex_type, int> const& combination)> const& func
    )
{
    if(it == conbination.end())
    {
        func(combination);
    }
    else
    {
        for(int i = 0; i < it->first->selectable_num; ++i)
        {
            it->second = i;
            all_combination_pattern(combination, it + 1, func);
        }
    }
}

// combinationから必要なノードの選択状態だけ取り出して，条件を更新する
std::unordered_map<vertex_type, int> update_select_condition(
    std::unordered_map<vertex_type, int> const& whole_condition,
    std::unordered_map<vertex_type, int> const& base_condition
    )
{
    std::unordered_map<vertex_type, int> particle_condition = base_condition;
    for(auto it = particle_condition.begin(); it != particle_condition.end(); ++it)
    {
        it->second = whole_condition[it->first];
    }
    return particle_condition;
}

// 周辺化
std::unordered_map<std::unordered_map<vertex_type, int>, double> marginalize(
    std::unordered_map<std::unordered_map<vertex_type, int>, double> const& base,
    vertex_type const& target
    )
{
    std::unordered_map<std::unordered_map<vertex_type, int>, double> result;
    all_combination_pattern(
        base, base.begin(),
        [&result, &base, &target](std::unordered_map<std::unordered_map<vertex_type, int>, double> const& condition)
        {
            auto new_condition = condition;
            new_condtion.erase(target);

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

likelihood_list lilelihood_list_;

// 条件化(条件を添加する)
void conditioning(
    std::unordered_map<std::unordered_map<vertex_type, int>, double> const& probabilities,
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

    for(int i = 0; i < condition_node.selectable_num; ++i)
    {
        double denominator = 0;
        for(int j = 0; j < target_node.selectable; ++j)
        {
            auto const target_prob = pair_probabilities[{{condition_node,i},{target_node,j}}];
            denominator += target_prob;
            lilelihood_list_(condition_node, target_node)[i][j] = target_prob;
        }

        for(int j = 0; j < target_node.selectable; ++j)
        {
            lilelihood_list_(condition_node, target_node)[i][j] /= denominator;
        }
    }
}

// 指定ノードから上流に拡散，上流が確定した後に自身を算出することで解を得る
std::unordered_map<std::unordered_map<vertex_type, int>, double> calculate_likelihood_from_backward(
    graph_t const& graph,
    vertex_type const& node
    )
{
    std::vector<vertex_type> direct_parents; // 直接の親
    std::vector<vertex_type> indirect_parents; // 間接の親

    std::vector<std::unordered_map<std::unordered_map<vertex_type, int>, double>> upward_probabilities;

    // 上流ノードの列挙(boost::transformを使うと簡潔)
    for(auto const& upward_edge : graph.in_edges(node))
    {
        auto const upward_node = graph.source(upward_edge);
        direct_parents.push_back(upward_node);

        likelihood_list_.ad_manage(upward_node, node, upward_edge);
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
            auto const& parent_node = probability->first;
            if(std::find(direct_parents.cbegin(), direct_parents.cend(), parent_node) == direct_parents.cend() &&
               std::find(indirect_parents.cbegin(), indirect_parents.cend(), parent_node) == indirect_parents.cend())
            {
                indirect_parents.push_back(parent_node);
            }
        }

        upward_probabilities.push_back(std::move(result));
    }

    // 上流ノード全ての組み合わせを作製して回す
    std::unordered_map<vertex_type, int> combination = {{node, 0}};
    for(auto const& key : direct_parents) combination[key] = 0;
    for(auto const& key : indirect_parents) combination[key] = 0;

    // returnにも使われる同時確率を計算
    std::unordered_map<std::unordered_map<vertex_type, int>, double> target_node_probability;
    all_combination_pattern(
        combination, combination.begin(),
        [&target_node_probability, &upward_probabilities](std::unordered_map<vertex_type, int> const& combination)
        {
            // foldl使うと，どうせVC落ちるからやめた
            double probability = node->cpt[update_select_condition(combination, node->cpt.begin()->first)]
                                          [combination[node]];

            for(auto const& upward : upward_probabilities)
            {
               probability *= upward[update_select_condition(combination, upward.begin()->first)];
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
void calculate_likelihood(graph_t const& graph)
{
    auto node_list = graph.vertex_list();
    for(auto const& node : node_list)
    {
        auto edges = graph.out_edges(node);
        if(edges.size() == 0)
        {
            // 末端から走査したほうが都合がいいので，末端のノードを全網羅(==全エッジ走査)
            calculate_likelihood_from_backward // TODO:
        }
    }
}

matrix_type bp::operator()(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
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
            if(!edge->likelihood.first) throw std::runtime_error("no set edge of likelihood");
            mat = mat % (edge->likelihood.second * propagate_forward(graph, graph.target(edge), condition));
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
            if(!edge->likelihood.first) throw std::runtime_error("no set edge of likelihood");
            mat = mat % (propagate_backward(graph, graph.source(edge), condition) * edge->likelihood.second);
        }

        return mat;
    }

    // 最上位ノードは事前確率を割り当てる
    auto& e = node->evidence;
    if(e.first)
    {
        return e.second;
    }
    else
    {
        throw std::runtime_error("highest node doesn't have prior probability.");
    }
}

} // namespace bn

