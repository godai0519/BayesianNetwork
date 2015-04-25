#ifndef BNI_LEARNING_GREEDY_HPP
#define BNI_LEARNING_GREEDY_HPP

#include <string>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/learning/utility.hpp>

namespace bn {
namespace learning {

template<class Eval>
class greedy {
public:
    greedy(std::string const& filepath)
        : filepath_(filepath), sampling_(filepath_), eval_(sampling_), engine_(make_engine<std::mt19937>())
    {
    }
    
    double operator()(graph_t& graph)
    {
        // graphの初期化
        graph.erase_all_edge();

        auto vertexes = graph.vertex_list();
        std::shuffle(vertexes.begin(), vertexes.end(), engine_);

        // bestな構造を保持しておく
        sampling_.load_sample(graph.vertex_list());
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

private:
    std::string filepath_;
    sampler sampling_;
    Eval eval_;
    std::mt19937 engine_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_GREEDY_HPP
