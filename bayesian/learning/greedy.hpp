#ifndef BNI_LEARNING_GREEDY_HPP
#define BNI_LEARNING_GREEDY_HPP

#include <string>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/utility.hpp>

namespace bn {
namespace learning {

template<class Eval>
class greedy {
public:
    greedy(bn::sampler const& sampling)
        : sampling_(sampling), eval_(sampling_), engine_(make_engine<std::mt19937>())
    {
    }

    double operator()(graph_t& graph)
    {
        return (*this)(graph, graph.vertex_list());
    }

    double operator()(graph_t& graph, std::vector<vertex_type> vertexes)
    {
        std::shuffle(vertexes.begin(), vertexes.end(), engine_);

        // bestな構造を保持しておく
        sampling_.make_cpt(graph);
        double eval_now = eval_(graph);

        for(auto it = vertexes.begin(); it != vertexes.end();)
        {
            auto const child_iter = it;
            std::shuffle(++it, vertexes.end(), engine_);

            for(auto parent_iter = it; parent_iter != vertexes.end(); ++parent_iter)
            {
                if(auto edge = graph.add_edge(*parent_iter, *child_iter))
                {
                    // 辺を貼れたなら調べてみる
                    sampling_.make_cpt(graph);
                    auto const eval_next = eval_(graph);

                    if(eval_next < eval_now)
                    {
                        // 良くなってる
                        eval_now = eval_next;
                    }
                    else
                    {
                        // 変わらない，もしくは悪い == 元に戻す
                        graph.erase_edge(edge);
                    }
                }
            }
        }

        return eval_now;
    }

    // ヒントを与えた上でgreedy探索を行う
    // graph: 元となったグラフ構造
    // parent_nodes, child_nodes: parent_nodesに含まれるnodeからchild_nodeに含まれるnodeにしか辺を張らない
    double learn_with_hint(graph_t& graph, std::vector<vertex_type> parent_nodes, std::vector<vertex_type> child_nodes)
    {
        // 子ノードの出現順序をランダムに
        std::shuffle(std::begin(child_nodes), std::end(child_nodes), engine_);

        // 最高評価値を保存しておく
        sampling_.make_cpt(graph);
        double eval_now = eval_(graph);

        // 親候補と子候補を全部回して様子見る
        for(auto const& child : child_nodes)
        {
            // 親ノードの出現順序をランダムに
            std::shuffle(std::begin(parent_nodes), std::end(parent_nodes), engine_);

            for(auto const& parent : parent_nodes)
            {
                if(auto edge = graph.add_edge(parent, child))
                {
                    // 辺が張れたならば，評価をする
                    sampling_.make_cpt(graph);
                    auto const eval_next = eval_(graph);

                    if(eval_next < eval_now)
                        // 良くなってるscore
                        eval_now = eval_next;
                    else
                        // 変わらない or 悪い -> 元に戻す
                        graph.erase_edge(edge);
                }
            }
        }

        return eval_now;
    }

private:
    sampler const& sampling_;
    Eval const eval_;
    std::mt19937 engine_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_GREEDY_HPP
