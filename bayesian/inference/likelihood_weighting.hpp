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
    typedef std::unordered_map<bn::condition_t, std::size_t> sample_list;
    typedef std::unordered_map<vertex_type, matrix_type> return_type;

    explicit likelihood_weighting(graph_t const& graph)
        : graph_(graph)
    {
    }

    virtual ~likelihood_weighting() = default;


    // Probability Inference with Likelihood Weighting. P(Q|E)
    // evidence:  evidence node(E) list
    // unit_size: {n * unit_size} samples are used (n = 1,2,...)
    // epsilon:   an acceptable error range
    return_type operator() (
        evidence_list const& evidence,
        std::size_t const unit_size = 1000000/* 1'000'000 */,
        double const epsilon = 0.001
        )
    {
        auto const& node_list = graph_.vertex_list();

        // Initialize
        return_type ret;
        for(auto const& node : graph_.vertex_list())
            ret[node].resize(1, node->selectable_num, 1.0 / node->selectable_num);

        // Patterns per Unit
        while(true)
        {
            // Copy Previous Unit
            return_type next = ret;

            // Generate Unit
            for(std::size_t i = 0; i < unit_size; ++i)
            {
                auto const sample = weighted_sample(evidence);
                auto const& pattern = sample.first;
                auto const& w = sample.second;

                for(auto const& node : node_list)
                {
                    auto const node_select = pattern.at(node);
                    next[node][0][node_select] += w;
                }
            }

            // Normalization and Max Difference
            double difference = std::numeric_limits<double>::min();
            for(auto const& node : node_list)
            {
                auto& now_probabilities = ret[node];
                auto& next_probabilities = next[node];
                next_probabilities = this->normalize(next_probabilities);

                for(std::size_t i = 0; i < node->selectable_num; ++i)
                {
                    difference = std::max(
                        difference,
                        std::abs(next_probabilities[0][i] - now_probabilities[0][i])
                        );
                }
            }

            // Next Unit ?
            ret = std::move(next);
            if(difference <= epsilon) break;
        }

        return ret;
    }

    // Make samples with Likelihood Weighting. P(Q|E)
    // evidence:   evidence node(E) list
    // sample_num: {unit_size} samples are used
    sample_list make_samples(
        evidence_list const& evidence,
        std::size_t const sample_num = 1000000/* 1'000'000 */
        )
    {
        sample_list samples;

        for(std::size_t i = 0; i < sample_num; ++i)
        {
            // Generate One
            auto const& sample = weighted_sample(evidence).first;

            // Find and Count up
            auto const it = samples.lower_bound(sample);
            if(it == samples.end() || it->first != sample)
            {
                samples.emplace_hint(
                    it,
                    std::piecewise_construct,
                    std::forward_as_tuple(sample),
                    std::forward_as_tuple(1)
                    );
            }
            else
            {
                ++(it->second);
            }
        }

        return samples;
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

        if(sum < 1.0e-20)
        {
            // 正規化出来ないレベルならば，一様にする
            auto const elem_num = target.height() * target.width();
            for(std::size_t i = 0; i < target.height(); ++i)
                for(std::size_t j = 0; j < target.width(); ++j)
                    target[i][j] = 1.00 / elem_num;
        }
        else
        {
            // 正規化できるならば，正規化する
            for(std::size_t i = 0; i < target.height(); ++i)
                for(std::size_t j = 0; j < target.width(); ++j)
                    target[i][j] /= sum;
        }
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

