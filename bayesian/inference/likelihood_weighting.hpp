#ifndef BNI_INFERENCE_LIKELIHOOD_WEIGHTING_HPP
#define BNI_INFERENCE_LIKELIHOOD_WEIGHTING_HPP

#include <algorithm>
#include <functional>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/matrix.hpp>

namespace bn {
namespace inference {

class likelihood_weighting {
public:
    typedef std::unordered_map<vertex_type, int> evidence_list;
    typedef std::unordered_map<vertex_type, int> pattern_list;
    typedef std::unordered_map<bn::condition_t, std::size_t> newer_pattern_list;
    typedef std::unordered_map<vertex_type, matrix_type> return_type;

    explicit likelihood_weighting(graph_t const& graph)
        : graph_(graph)
    {
    }

    virtual ~likelihood_weighting() = default;

    // Run: Likelihood Weighting
    return_type operator() (evidence_list const& evidence, std::uint64_t const sample_num = 10000)
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
            ret[node] = normalize(ret[node]);
        }

        return ret;
    }

    // Make: accurate sample
    std::pair<newer_pattern_list, return_type> make_samples_min(
        evidence_list const& evidence,
        std::uint64_t const unit_size/* = 1000000/* 1'000'000 */,
        double const epsilon/* = 0.001*/
        )
    {
        newer_pattern_list patterns;

        return_type w_list;
        return_type probabilities;

        // Initialize
        for(auto const& node : graph_.vertex_list())
        {
            w_list[node].resize(1, node->selectable_num, 0.0);
            probabilities[node].resize(1, node->selectable_num, 0.0);
        }

        while(true)
        {
            // Generate one unit
            for(std::uint64_t i = 0; i < unit_size; ++i)
            {
                auto const sample = weighted_sample(evidence);
                pattern_list const& pattern = sample.first;
                double const& w = sample.second;

                for(auto const& node : graph_.vertex_list())
                {
                    auto const node_select = pattern.at(node);
                    w_list[node][0][node_select] += w;
                }

                auto it = patterns.find(pattern);
                if(it != patterns.cend()) ++(it->second);
                else                      patterns[pattern] = 1;
            }

            double max_difference = std::numeric_limits<double>::min();
            for(auto const& node : graph_.vertex_list())
            {
                auto& old_probability = probabilities[node];
                auto next_probability = normalize(w_list[node]);
                for(std::size_t i = 0; i < next_probability.width(); ++i)
                {
                    max_difference = std::max(max_difference, std::abs(old_probability[0][i] - next_probability[0][i]));
                }

                old_probability = std::move(next_probability);
            }

            if(max_difference < epsilon) break;
        }

        return std::make_pair(patterns, probabilities);
    }

    // Make: accurate sample
    // 非推奨化する
    std::pair<std::vector<pattern_list>, return_type> make_samples(
        evidence_list const& evidence,
        std::uint64_t const unit_size/* = 1000000/* 1'000'000 */,
        double const epsilon/* = 0.001*/,
        std::function<void(std::vector<pattern_list>&)> handler
        )
    {
        std::vector<pattern_list> patterns;
        return_type w_list;
        return_type probabilities;

        // Initialize
        for(auto const& node : graph_.vertex_list())
        {
            w_list[node].resize(1, node->selectable_num, 0.0);
            probabilities[node].resize(1, node->selectable_num, 0.0);
        }

        while(true)
        {
            // Generate one unit
            patterns.reserve(unit_size + patterns.capacity());
            for(std::uint64_t i = 0; i < unit_size; ++i)
            {
                auto const sample = weighted_sample(evidence);
                pattern_list const& pattern = sample.first;
                double const& w = sample.second;

                for(auto const& node : graph_.vertex_list())
                {
                    auto const node_select = pattern.at(node);
                    w_list[node][0][node_select] += w;
                }

                patterns.push_back(pattern);
            }

            handler(patterns);

            double max_difference = std::numeric_limits<double>::min();
            for(auto const& node : graph_.vertex_list())
            {
                auto& old_probability = probabilities[node];
                auto next_probability = normalize(w_list[node]);
                for(std::size_t i = 0; i < next_probability.width(); ++i)
                {
                    max_difference = std::max(max_difference, std::abs(old_probability[0][i] - next_probability[0][i]));
                }

                old_probability = std::move(next_probability);
            }

            if(max_difference < epsilon) break;
        }

        return std::make_pair(patterns, probabilities);
    }

private:
    // Make a sample according to evidence node
    // evidenceノードに従って，1サンプルを作成する
    std::pair<pattern_list, double> weighted_sample(evidence_list const& evidence)
    {
        double w = 1.0;
        pattern_list pattern;

        // 再帰関数
        std::function<void(graph_t const&, vertex_type const&, std::vector<vertex_type>&)> recursion;
        recursion = [this, &recursion, &w, &pattern, &evidence](graph_t const& graph_, vertex_type const& target, std::vector<vertex_type>& remain_node)
            {
                condition_t parent_cond;

                for(auto const& parent : graph_.in_vertexes(target))
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

    // Using weight, select from random value
    // weightを使って，ランダム値からそれがどの選択値になるか判断する
    int make_random_by_weight(double const value, std::vector<double> const& weight) const
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

    // Normalization Σ[i,j](a_ij) = 1
    // (重複)
    matrix_type normalize(matrix_type target) const
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

    // 乱数生成器
    class probability_generator {
    public:
        probability_generator()
            : distribution_(0.0, 1.0)
        {
            std::random_device rand_dev;
            std::vector<std::uint_least32_t> vec(10);
            std::generate(vec.begin(), vec.end(), std::ref(rand_dev));
            std::seed_seq seed(vec.begin(), vec.end());
            engine_.reset(new std::mt19937(seed));
        }

        double operator() ()
        {
            return distribution_(*engine_);
        }

    private:
        std::unique_ptr<std::mt19937> engine_;
        std::uniform_real_distribution<double> distribution_;
    };

    graph_t const graph_;
    probability_generator probability_generator_;
};

} // namespace inference
} // namespace bn

#endif // #ifndef BNI_INFERENCE_LIKELIHOOD_WEIGHTING_HPP

