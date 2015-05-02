#ifndef BNI_LEARNING_BRUTE_FORCE_HPP
#define BNI_LEARNING_BRUTE_FORCE_HPP

#include <string>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/learning/utility.hpp>

namespace bn {
namespace learning {

template<class Eval>
class brute_force {
public:
    brute_force(bn::sampler const& sampling)
        : sampling_(sampling), eval_(sampling_)
    {
    }

    // グラフに属するすべての頂点について，全探索を行い，最適解にする
    // 最適解をグラフ構造，最適スコアをreturnする
    double operator()(graph_t& graph)
    {
        return (*this)(graph, graph.vertex_list());
    }

    // 指定した頂点について，全探索を行う
    // 最適解をグラフ構造，最適スコアをreturnする
    // graph: グラフ構造
    // vertexes: 学習する頂点
    double operator()(graph_t& graph, std::vector<vertex_type> const& vertexes)
    {
        // bestな構造を保持しておく
        sampling_.make_cpt(graph);
        graph_t best_graph = graph;
        double best_eval  = eval_(graph);

        // 全探索する
        recursive(0, graph, vertexes, best_graph, best_eval);

        graph = std::move(best_graph);
        return best_eval;
    }

    // 指定した親ノード群・子ノード群に対し，この条件に沿った構造を全探索する
    // 最適解をグラフ構造，最適スコアをreturnする
    // graph: グラフ構造
    // parent_nodes: 親ノードにしかなれないノードの集合
    // child_nodes: 子ノードにしかなれないノードの集合
    double learn_with_hint(graph_t& graph, std::vector<vertex_type> parent_nodes, std::vector<vertex_type> child_nodes)
    {
        // 可能な辺の数
        std::size_t const possible_edge_num = parent_nodes.size() * child_nodes.size();

        // 可能な辺のリスト
        std::vector<std::pair<vertex_type, vertex_type>> possible_edge;
        possible_edge.reserve(possible_edge_num);

        // リストアップ
        for(auto const& parent : parent_nodes)
        {
            for(auto const& child : child_nodes)
            {
                possible_edge.emplace_back(parent, child);
            }
        }

        // bestな構造を保持しておく
        sampling_.make_cpt(graph);
        graph_t best_graph = graph;
        double best_eval  = eval_(graph);

        // 全辺の全通りについて探索
        recursive_with_hint(graph, possible_edge.cbegin(), possible_edge.cend(), best_graph, best_eval);

        // 最終処理
        graph = std::move(best_graph);
        return best_eval;
    }

private:
    // learn_with_hint用再帰関数
    template<class Iterator>
    void recursive_with_hint(
        graph_t& graph, Iterator const& begin, Iterator const& end,
        graph_t& best_graph, double& best_eval
        )
    {
        if(begin == end)
        {
            // 評価
            sampling_.make_cpt(graph);
            auto const now_eval = eval_(graph);
            if(now_eval < best_eval)
            {
                best_eval = now_eval;
                best_graph = graph;
            }
        }
        else
        {
            // 張らない
            recursive_with_hint(graph, begin + 1, end, best_graph, best_eval);

            // 張る
            if(auto const edge = graph.add_edge(begin->first, begin->second))
            {
                recursive_with_hint(graph, begin + 1, end, best_graph, best_eval);
                graph.erase_edge(edge);
            }
        }
    }

    // operator()用再帰関数
    void recursive(
        std::size_t const& target,
        graph_t& graph, std::vector<vertex_type> const& vertexes,
        graph_t& best_graph, double& best_eval
        )
    {
        // 情報整理
        auto const vertex_num = vertexes.size();

        if(target == vertex_num - 1)
        {
            // 評価
            sampling_.make_cpt(graph);
            auto const now_eval = eval_(graph);
            if(now_eval < best_eval)
            {
                best_eval = now_eval;
                best_graph = graph;
            }
        }
        else
        {
            // target <-> i 間に，「張らない」「右向き」「左向き」
            for(std::size_t i = target + 1; i < vertex_num; ++i)
            {
                recursive(target + 1, graph, vertexes, best_graph, best_eval);

                if(auto const edge = graph.add_edge(vertexes[target], vertexes[i]))
                {
                    recursive(target + 1, graph, vertexes, best_graph, best_eval);
                    graph.erase_edge(edge);
                }

                if(auto const edge = graph.add_edge(vertexes[i], vertexes[target]))
                {
                    recursive(target + 1, graph, vertexes, best_graph, best_eval);
                    graph.erase_edge(edge);
                }
            }
        }
    }

    sampler const& sampling_;
    Eval const eval_;
};


} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_BRUTE_FORCE_HPP
