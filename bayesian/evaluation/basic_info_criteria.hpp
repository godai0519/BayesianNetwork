#ifndef BNI_EVALUATION_BASIC_INFO_CRITERIA_HPP
#define BNI_EVALUATION_BASIC_INFO_CRITERIA_HPP

#include <cstdint>
#include <numeric>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/evaluation/basic_evaluation.hpp>

namespace bn {
namespace evaluation {

class basic_info_criteria : basic_evaluation{
public:
    basic_info_criteria(sampler const& sampling);

    // API
    virtual double operator() (graph_t const& graph) const
    {
        return (*this)(graph, graph.vertex_list());
    }
    virtual double operator() (graph_t const& graph, std::vector<bn::vertex_type> const& vertex_list) const = 0;

protected:
    // AIC/MDL情報基準共通: ネットワーク構造の適切さに関する項
    // - log P_theta^N(D)
    double calc_likelihood(graph_t const& graph) const;
    double calc_likelihood(graph_t const& graph, std::vector<bn::vertex_type> const& vertex_list) const;

    // AIC/MDL情報基準共通: ネットワーク構造の複雑さに関する項
    // + d
    double calc_parameters(graph_t const& graph) const;

    // Sampling数のgetter
    std::size_t sampling_size() const;

private:
    sampler const& sampling_;
};

basic_info_criteria::basic_info_criteria(sampler const& sampling)
    : sampling_(sampling)
{
}

double basic_info_criteria::calc_likelihood(graph_t const& graph) const
{
    return calc_likelihood(graph, graph.vertex_list());
}

double basic_info_criteria::calc_likelihood(graph_t const& graph, std::vector<bn::vertex_type> const& vertex_list) const
{
    double likelihood = 0.0;

    // collect sample data to statistics
    std::unordered_map<vertex_type, std::unordered_map<condition_t, std::vector<std::size_t>>> statistics;
    for(auto const& sample : sampling_.table())
    {
        for(auto const& node : vertex_list)
        {
            condition_t cond;
            for(auto const& parent : graph.in_vertexes(node))
                cond[parent] = sample.first.at(parent);

            auto& table = statistics[node];
            auto it = table.lower_bound(cond);
            if(it == table.end() || it->first != cond)
            {
                it = table.emplace_hint(
                    it, std::piecewise_construct,
                    std::forward_as_tuple(cond),
                    std::forward_as_tuple(node->selectable_num, 0)
                    );
            }

            it->second[sample.first.at(node)] += sample.second;
        }
    }

    // from statistics, calculating log-likelihood 
    for(auto const& node_statistics : statistics)
    {
        for(auto it = node_statistics.second.begin(); it != node_statistics.second.end(); ++it)
        {
            std::size_t const sum = std::accumulate(it->second.begin(), it->second.end(), 0);
            for(std::size_t i = 0; i < node_statistics.first->selectable_num; ++i)
            {
                auto const score = std::log2(static_cast<double>(it->second[i]) / sum);
                if(!std::isinf(score))
                    likelihood -= it->second[i] * score;
            }
        }
    }

    return likelihood;
}

double basic_info_criteria::calc_parameters(graph_t const& graph) const
{
    int64_t parameters = 0;
    for(auto const& node : graph.vertex_list())
    {
        auto const& parents = graph.in_vertexes(node);

        // 1ノードのパラメータ数を足し合わせる
        parameters += std::accumulate(
            parents.begin(), parents.end(), node->selectable_num - 1,
            [](int64_t const& init, bn::vertex_type const& parent_node) -> int64_t
            {
                return init * parent_node->selectable_num;
            });
    }

    return static_cast<double>(parameters);
}

std::size_t basic_info_criteria::sampling_size() const
{
    return sampling_.sampling_size();
}

} // namespace evaluation
} // namespace bn

#endif // #ifndef BNI_EVALUATION_BASIC_INFO_CRITERIA_HPP
