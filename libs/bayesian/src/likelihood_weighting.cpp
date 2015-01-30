#include <algorithm>
#include <functional>
#include "bayesian/graph.hpp"
#include "bayesian/likelihood_weighting.hpp"

namespace bn {

likelihood_weighting::likelihood_weighting(graph_t const& graph)
    : graph_(graph)
{
}

auto likelihood_weighting::operator() (evidence_list const& evidence, std::uint64_t const sample_num) -> return_type
{
    return_type ret;

    // Initialize
    for(auto const& node : graph_.vertex_list())
    {
        ret[node].resize(1, node->selectable_num, 0.0);
    }

    // Generate Sample
    for(std::uint64_t i = 0; i < sample_num; ++i)
    {
        auto const sample = weighted_sample(evidence);
        pattern_list const& pattern = sample.first;
        double const& w = sample.second;
        
        for(auto const& node : graph_.vertex_list())
        {
            auto const node_select = pattern.at(node);
            ret[node][0][node_select] += w;
        }
    }
    
    // Normalization
    for(auto const& node : graph_.vertex_list())
    {
        normalize(ret[node]);
    }

    return ret;
}

auto likelihood_weighting::weighted_sample(evidence_list const& evidence) -> std::pair<pattern_list, double>
{
    double w = 1.0;
    pattern_list pattern;

    // 再帰関数
    std::function<void(graph_t const&, vertex_type const&, std::vector<vertex_type>&)> recursion;
    recursion = [this, &recursion, &w, &pattern, &evidence](graph_t const& graph_, vertex_type const& target, std::vector<vertex_type>& remain_node)
        {
            condition_t parent_cond;

            for(auto const& parent : graph_.in_vertexs(target))
            {
                // 親がまだ設定されていなければ再帰
                auto const it = std::find(remain_node.cbegin(), remain_node.cend(), parent);
                if(it != remain_node.cend())
                {
                    remain_node.erase(it);
                    recursion(graph_, parent, remain_node);
                }

                // 条件追加
                parent_cond[parent] = pattern[parent];
            }
            
            // evidenceに含まれれば or 含まれなければ乱数
            auto const evidence_it = evidence.find(target);
            if(evidence_it != evidence.end())
            {
                w *= target->cpt[parent_cond].second.at(evidence_it->second);
                pattern[target] = evidence_it->second;
            }
            else
            {
                int const select = make_random_by_weight(probability_generator_(), target->cpt[parent_cond].second);
                pattern[target] = select;
            }
        };

    // ノードがなくなるまで再帰関数を回す
    auto node_list = graph_.vertex_list();
    while(!node_list.empty())
    {
        // 最後尾のノードを対象とする
        auto const node = node_list.back();
        node_list.pop_back();

        recursion(graph_, node, node_list);
    }

    return std::make_pair(pattern, w);
}

int likelihood_weighting::make_random_by_weight(double const value, std::vector<double> const& weight) const
{
    assert(0.0 <= value && value < 1.0);

    double total = 0.0;
    for(int i = 0; i < weight.size(); ++i)
    {
        auto const old_total = total;
        total += weight.at(i);
        if(old_total <= value && value < total)
        {
            return i;
        }
    }

    return weight.size() - 1;
}

matrix_type& likelihood_weighting::normalize(matrix_type& target) const
{
    double sum = 0;

    for(std::size_t i = 0; i < target.height(); ++i)
        for(std::size_t j = 0; j < target.width(); ++j)
            sum += target[i][j];

    for(std::size_t i = 0; i < target.height(); ++i)
        for(std::size_t j = 0; j < target.width(); ++j)
            target[i][j] /= sum;

    return target;
}

likelihood_weighting::probability_generator::probability_generator()
    : distribution_(0.0, 1.0)
{
    std::random_device rand_dev;
    std::vector<std::uint_least32_t> vec(10);
    std::generate(vec.begin(), vec.end(), std::ref(rand_dev));
    std::seed_seq seed(vec.begin(), vec.end());
    engine_.reset(new std::mt19937(seed));
}

double likelihood_weighting::probability_generator::operator() ()
{
    return distribution_(*engine_);
}

} // namespace bn

