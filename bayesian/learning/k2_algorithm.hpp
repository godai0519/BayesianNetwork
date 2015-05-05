#ifndef BNI_LEARNING_K2_ALGORITHM_HPP
#define BNI_LEARNING_K2_ALGORITHM_HPP

#include <string>
#include <random>
#include <boost/range/algorithm_ext/erase.hpp>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/learning/utility.hpp>

namespace bn {
namespace learning {
    
template<class Eval>
class k2_algorithm {
public:
    k2_algorithm(bn::sampler const& sampling)
        : sampling_(sampling), eval_(sampling_), engine_(make_engine<std::mt19937>())
    {
    }
    
    double operator()(graph_t& graph, std::unordered_map<vertex_type, std::vector<vertex_type>> preconditon)
    {   
        // bestな構造を保持しておく
        sampling_.make_cpt(graph);
        double eval_best = eval_(graph);

        // 頂点の順番をランダムに
        auto vertexes = graph.vertex_list();
        std::shuffle(vertexes.begin(), vertexes.end(), engine_);

        for(auto const& target : vertexes)
        {
            // 自分自身もしくは与えられた除外条件以外のノードを親候補とする
            auto candidature = graph.vertex_list();
            boost::remove_erase_if(candidature,
                [&target, &preconditon](vertex_type const& node) -> bool
                {
                    // C言語でいう副作用完了点を悪用
                    auto ignore_nodes = preconditon.find(target);
                    return 
                        node == target || 
                        (ignore_nodes != preconditon.end() &&
                        std::find(ignore_nodes->second.begin(), ignore_nodes->second.end(), node) != ignore_nodes->second.end());
                });

            for(auto const& parent : candidature)
            {
                if(auto edge = graph.add_edge(parent, target))
                {
                    sampling_.make_cpt(graph);
                    auto const eval_now = eval_(graph);
                    
                    if(eval_now < eval_best)
                    {
                        eval_best = eval_now;
                        preconditon[parent].push_back(target);
                    }
                    else graph.erase_edge(edge);
                }
            }
        }

        return eval_best;
    }

private:
    sampler const& sampling_;
    Eval const eval_;
    std::mt19937 engine_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_K2_ALGORITHM_HPP
