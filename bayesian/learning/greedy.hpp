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
        // CPT生成
        sampling_.make_cpt(graph);

        // 子候補とし全て回す
        std::shuffle(vertexes.begin(), vertexes.end(), engine_); // 子ノードの出現順序をランダムに
        for(auto it = vertexes.begin(); it != vertexes.end();)
        {
            // 子ノードの評価値を得る（親ノードは評価値は変化しない）
            auto const child_iter = it;
            double evaluation_value = eval_(graph, {*child_iter});

            // 親ノードすべてに結合テストを行う
            std::shuffle(++it, vertexes.end(), engine_); // 親ノードの出現順序をランダムに
            for(auto parent_iter = it; parent_iter != vertexes.end(); ++parent_iter)
            {
                if(auto edge = graph.add_edge(*parent_iter, *child_iter))
                {
                    // 辺が張れたならば，評価をする
                    auto old_cpt = (*child_iter)->cpt; // Backup
                    sampling_.make_cpt(graph, *child_iter);
                    double const next_evaluation_value = eval_(graph, {*child_iter});
                    
                    if(next_evaluation_value < evaluation_value)
                    {
                        // 良くなってる
                        evaluation_value = next_evaluation_value;
                    }
                    else
                    {
                        // 変わらない，もしくは悪い == 元に戻す
                        graph.erase_edge(edge);
                        (*child_iter)->cpt = old_cpt; // Restore
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
        // CPT生成
        sampling_.make_cpt(graph);

        // 親候補と子候補を全部回して様子見る
        std::shuffle(std::begin(child_nodes), std::end(child_nodes), engine_); // 子ノードの出現順序をランダムに
        for(auto const& child : child_nodes)
        {
            // 子ノードの評価値を得る（親ノードは評価値は変化しない）
            double evaluation_value = eval_(graph, {child});
            
            // 親ノードすべてに結合テストを行う
            std::shuffle(std::begin(parent_nodes), std::end(parent_nodes), engine_); // 親ノードの出現順序をランダムに
            for(auto const& parent : parent_nodes)
            {
                if(auto edge = graph.add_edge(parent, child))
                {
                    // 辺が張れたならば，評価をする
                    auto old_cpt = child->cpt; // Backup
                    sampling_.make_cpt(graph, child);
                    double const next_evaluation_value = eval_(graph, {child});

                    if(next_evaluation_value < evaluation_value)
                    {
                        // 良くなってるscore
                        evaluation_value = next_evaluation_value;
                    }
                    else
                    {
                        // 変わらない or 悪い -> 元に戻す
                        graph.erase_edge(edge);
                        child->cpt = old_cpt; // Restore
                    }
                }
            }
        }

        return std::numeric_limits<double>::quiet_NaN();
    }

private:
    sampler const& sampling_;
    Eval const eval_;
    std::mt19937 engine_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_GREEDY_HPP
