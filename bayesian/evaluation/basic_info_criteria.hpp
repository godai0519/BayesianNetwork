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
    virtual double operator() (graph_t const& graph) const = 0;

protected:
    // AIC/MDL情報基準共通: ネットワーク構造の適切さに関する項
    // - log P_theta^N(D)
    double calc_likelihood(graph_t const& graph) const;

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
    double likelihood = 0.0;

    for(auto const& sample : sampling_.table())
    {
        // 各ノードに対し
        for(auto const& node : graph.vertex_list())
        {
            condition_t cond = sample.first;

            // 親ノード以外を消す
            auto const& parent = graph.in_vertexes(node);
			for (auto it = cond.begin(); it != cond.end();)
			{
				if(std::find(parent.begin(), parent.end(), it->first) == parent.end())
				{
					it = cond.erase(it);
				}
				else
				{
					++it;
				}
			}

            // 全選択について計算
			for (int i = 0; i < node->selectable_num; i++)
            {
				if (node->cpt[cond].second[i] != 0)
                {
					likelihood -= sample.second * std::log2(node->cpt[cond].second[i]);
				}
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
                return init * (parent_node->selectable_num - 1);
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
