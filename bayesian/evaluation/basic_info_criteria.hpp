#ifndef BNI_EVALUATION_BASIC_INFO_CRITERIA_HPP
#define BNI_EVALUATION_BASIC_INFO_CRITERIA_HPP

#include <cstdint>
#include <numeric>
#include <boost/optional.hpp>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/evaluation/basic_evaluation.hpp>

namespace bn {
namespace evaluation {

class basic_info_criteria : basic_evaluation{
public:
    basic_info_criteria(std::string const& file);

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
    boost::optional<std::size_t> const& sampling_size() const;

private:
    sampler const loader_;
    std::string const file_;
};

basic_info_criteria::basic_info_criteria(std::string const& file)
    : loader_(), file_(file)
{
}

double basic_info_criteria::calc_likelihood(graph_t const& graph) const
{
    double likelihood = 0.0;

    // サンプルを1行ずつ読み込む
    loader_.load_sample(
        graph, file_,
        [&likelihood, &graph](condition_t const& sample)
        {
            condition_t cond = sample;

            // 各ノードに対し
            for(auto const& node : graph.vertex_list())
            {
				for (auto it = cond.begin(); it != cond.end();)
				{
					auto const& parent = graph.in_vertexes(node);
					if (std::find(parent.begin(), parent.end(), it->first) == parent.end())
					{
						it = cond.erase(it);
					}
					else
					{
						++it;
					}
				}
				for (int i = 0; i < node->selectable_num; i++)
                {
					if (node->cpt[cond].second[i] != 0)
                    {
						likelihood -= std::log2(node->cpt[cond].second[i]);
					}
				}
            }
        });

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

boost::optional<std::size_t> const& basic_info_criteria::sampling_size() const
{
    return loader_.sampling_size();
}

} // namespace evaluation
} // namespace bn

#endif // #ifndef BNI_EVALUATION_BASIC_INFO_CRITERIA_HPP
