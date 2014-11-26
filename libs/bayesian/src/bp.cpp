#include <boost/foreach.hpp>
#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"

#include <iostream>

namespace bn {

matrix_type bp::operator()(
    graph_t const& graph,
    graph_t::vertex_descriptor const& node,
    std::vector<std::pair<graph_t::vertex_descriptor, int>> const& evidence
    )
{
    // 前後の要素に伝播させる
    auto const e_minus = propagate_forward(graph, node, evidence);
    auto const e_plus = propagate_backward(graph, node, evidence);

    // 掛け算
    auto const elem_num = graph[node].selectable_num;
    double sum = 0.0;
    matrix_type mat(boost::extents[elem_num][1]);
    for(int i = 0; i < e_minus.shape()[0]; ++i)
    {
        double const product = e_minus[i][0] * e_plus[0][i];
        sum += product;
        mat[i][0] = product;
    }

    // 正規化
    for(int i = 0; i < e_minus.shape()[0]; ++i)
    {
        mat[i][0] /= sum;
    }

    return mat;
}

std::pair<bool, int> find_evidence(
    graph_t::vertex_descriptor const& node,
    std::vector<std::pair<graph_t::vertex_descriptor, int>> const& evidence
    )
{
    for(auto const& e : evidence)
    {
        if(e.first == node)
        {
            return std::make_pair(true, e.second);
        }
    }

    return std::make_pair(false, 0);
}

// 下流要素の確率推論
matrix_type bp::propagate_forward(
    graph_t const& graph,
    graph_t::vertex_descriptor const& node,
    std::vector<std::pair<graph_t::vertex_descriptor, int>> const& evidence
    )
{
    // node ∈ evidence
    auto is_evidence =  find_evidence(node, evidence);
    if(is_evidence.first)
    {
        auto const elem_num = graph[node].selectable_num;
        matrix_type mat(boost::extents[elem_num][1]);
        for(int i = 0; i < elem_num; ++i)
        {
            mat[i][0] = (i == is_evidence.second) ? 1 : 0;
        }
        return mat;
    }

    // evidenceでないから伝播 (e-要素)
    BOOST_FOREACH(auto const& edge, out_edges(node, graph))
    {
        // TODO: 2つ以上の要素があった場合の処理を切り分ける
        return (*graph[edge].likelihood) * propagate_forward(graph, target(edge, graph), evidence);
    }

    // 末端は全ての確率が等しいとする
    auto const elem_num = graph[node].selectable_num;
    matrix_type mat(boost::extents[elem_num][1]);
    for(int i = 0; i < elem_num; ++i)
    {
        mat[i][0] = 1.0;
    }
    return mat;
}

// 上流要素の確率推論
matrix_type bp::propagate_backward(
    graph_t const& graph,
    graph_t::vertex_descriptor const& node,
    std::vector<std::pair<graph_t::vertex_descriptor, int>> const& evidence
    )
{
    // node ∈ evidence
    auto is_evidence =  find_evidence(node, evidence);
    if(is_evidence.first)
    {
        auto const elem_num = graph[node].selectable_num;
        matrix_type mat(boost::extents[1][elem_num]);
        for(int i = 0; i < elem_num; ++i)
        {
            mat[0][i] = (i == is_evidence.second) ? 1 : 0;
        }
        return mat;
    }

    // evidenceでないから伝播 (e+要素)
    BOOST_FOREACH(auto const& edge, in_edges(node, graph))
    {
        // TODO: 2つ以上の要素があった場合の処理を切り分ける
        return propagate_backward(graph, source(edge, graph), evidence) * (*graph[edge].likelihood);
    }

    // 末端は全てが等しいとする
    auto const elem_num = graph[node].selectable_num;
    matrix_type mat(boost::extents[1][elem_num]);
    for(int i = 0; i < elem_num; ++i)
    {
        mat[0][i] = 1.0;
    }
    return mat;
}

} // namespace bn

