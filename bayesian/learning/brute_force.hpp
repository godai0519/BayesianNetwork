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

    double operator()(graph_t& graph)
    {
        return (*this)(graph, graph.vertex_list());
    }

    double operator()(graph_t& graph, std::vector<vertex_type> const& vertexes)
    {
        // bestな構造を保持しておく
        sampling_.make_cpt(graph);
        graph_t best_graph = graph;
        double best_eval  = eval_(graph);;
        
        // 全探索する
        recursive(0, graph, vertexes, best_graph, best_eval);

        graph = std::move(best_graph);
        return best_eval;
    }

private:
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
