#ifndef BNI_LEARNING_SIMULATED_ANNEALING_HPP
#define BNI_LEARNING_SIMULATED_ANNEALING_HPP

#include <string>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/learning/utility.hpp>

namespace bn {
namespace learning {
    
template<class Eval>
class simulated_annealing {
public:
    simulated_annealing(bn::sampler const& sampling)
        : sampling_(sampling), eval_(sampling_), engine_(make_engine<std::mt19937>()), method_dist_(0, 2), probability_dist_(0.0, 1.0)
    {
    }

    double operator()(
        graph_t& graph,
        double const initial_temp, double const final_temp,
        double const decreasing_rate,
        double const boltzmann = 1.0,
        unsigned int const same_state_max = 100
        )
    {
        // bestな構造を保持しておく
        sampling_.make_cpt(graph);
        graph_t best_graph = graph;
        double best_eval = eval_(graph);

        // 情報整理
        auto const vertexes = graph.vertex_list();
        auto const vertex_num = vertexes.size();
        std::uniform_int<std::size_t> vertex_dist(0, vertex_num - 1);

        unsigned int no_changed_num = 0;
        double temperature = initial_temp;
        while(temperature >= final_temp && no_changed_num < same_state_max)
        {
            bool is_operated = false;
            auto const method = method_dist_(engine_);
            if(method == 0)
            {
                //
                // 辺の追加
                //

                auto const from = vertexes[vertex_dist(engine_)];
                auto const to   = vertexes[vertex_dist(engine_)];

                // Operate!
                if(graph.add_edge(from, to)) is_operated = true;
            }
            else if(method == 1)
            {
                //
                // 辺の削除
                //

                auto const edges = graph.edge_list();
                if(edges.size() < 1) continue;
                std::uniform_int<std::size_t> edge_dist(0, edges.size() - 1);

                // Operate!
                auto const target_edge = edges[edge_dist(engine_)];
                if(graph.erase_edge(target_edge)) is_operated = true;
            }
            else if(method == 2)
            {
                //
                // 辺の向きの変更
                //

                auto const edges = graph.edge_list();
                if(edges.size() < 1) continue;
                std::uniform_int<std::size_t> edge_dist(0, edges.size() - 1);

                // Operate!
                auto const target_edge = edges[edge_dist(engine_)];
                if(graph.change_edge_direction(target_edge)) is_operated = true;
            }

            // 何も変化していないのは，受理確率以前の問題なので，何かができるまで繰り返す
            if(!is_operated) continue;

            // 評価
            sampling_.make_cpt(graph);
            double const now_eval = eval_(graph);
            double const diff_eval = now_eval - best_eval;

            // 受理されたかどうか
            bool is_acceptance;
            if(diff_eval <= 0) is_acceptance = true;
            else               is_acceptance = probability_dist_(engine_) < std::exp(-now_eval / (boltzmann * temperature));

            if(is_acceptance)
            {
                // 移行
                best_graph = graph;
                best_eval = now_eval;
                no_changed_num = 0;
            }
            else
            {
                // 引き戻し
                graph = best_graph;
                ++no_changed_num;
            }

            // 焼きなまし
            temperature *= decreasing_rate;
        }

        return best_eval;
    }
    
private:
    sampler const& sampling_;
    Eval const eval_;
    std::mt19937 engine_;
    std::uniform_int_distribution<int> method_dist_;
    std::uniform_real_distribution<double> probability_dist_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_SIMULATED_ANNEALING_HPP
