#include <boost/foreach.hpp>
#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"

namespace bn {

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
    // node ∈ condition
    auto is_condition = find_condition(node, condition);
    if(is_condition.first)
    {
        auto const elem_num = node->selectable_num;
        matrix_type mat(elem_num, 1);
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
        auto const elem_num = node->selectable_num;
        matrix_type mat(elem_num, 1, 1);

        BOOST_FOREACH(auto const& edge, out_edges)
        {
            if(!edge->likelihood.first) throw std::runtime_error("no set edge of likelihood");
            mat = mat % (edge->likelihood.second * propagate_forward(graph, graph.target(edge), condition));
        }

        return mat;
    }

    // 末端は全ての確率が等しいとする
    auto const elem_num = node->selectable_num;
    matrix_type mat(elem_num, 1);
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
    // node ∈ condition
    auto is_condition = find_condition(node, condition);
    if(is_condition.first)
    {
        auto const elem_num = node->selectable_num;
        matrix_type mat(1, elem_num);
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
        auto const elem_num = node->selectable_num;
        matrix_type mat(1, elem_num, 1);

        BOOST_FOREACH(auto const& edge, in_edges)
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

