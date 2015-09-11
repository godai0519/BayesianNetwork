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
    struct element_type {
        // TODO:
        element_type() = default;
        element_type(std::vector<std::size_t> select, std::size_t num)
            : select(std::move(select)), num(num)
        {
        }

        std::vector<std::size_t> select;
        std::size_t              num;
    };

    typedef std::unordered_map<vertex_type, std::size_t> evidence_list;
    typedef std::vector<std::size_t>  pattern_list;
    typedef std::vector<element_type> sample_list;
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
        auto const topological = graph_.topological_sort();

        // Initialize
        return_type ret;
        for(auto const& node : graph_.vertex_list())
            ret[node].resize(1, node->selectable_num, 0.0);

        // Patterns per Unit
        while(true)
        {
            // Copy Previous Unit
            return_type next = ret;

            // Generate Unit
            for(std::size_t i = 0; i < unit_size; ++i)
            {
                auto const sample = weighted_sample(topological, evidence);
                auto const& pattern = sample.first;
                auto const& w = sample.second;

                for(auto const& node : node_list)
                {
                    auto const node_select = pattern.at(index_node_in_graph(node));
                    next[node][0][node_select] += w;
                }
            }

            // Normalization and Max Difference
            double difference = std::numeric_limits<double>::min();
            for(auto const& node : node_list)
            {
                auto const now_probabilities = normalize(ret[node]);
                auto const next_probabilities = normalize(next[node]);

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
        
        // Normalize
        for(auto const& node : node_list)
            ret[node] = normalize(ret[node]);

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
        auto const topological = graph_.topological_sort();

        for(std::size_t i = 0; i < sample_num; ++i)
        {
            // Generate One
            auto const& sample = weighted_sample(topological, evidence).first;

            // Find and Count up
            auto const it = std::find_if(
                samples.begin(), samples.end(),
                [&sample](element_type const& elem){ return elem.select == sample; }
                );

            if(it == samples.end())
                samples.emplace_back(sample, 1);
            else
                it->num += 1;
        }

        return samples;
    }

private:
    // Make a sample according to evidence node
    // evidenceノードに従って，1サンプルを作成する
    std::pair<pattern_list, double> weighted_sample(std::vector<vertex_type> const& topological, evidence_list const& evidence)
    {
        double w = 1.0;
        pattern_list pattern(topological.size());

        for(auto const& node : topological)
        {
            // Make CPT's condition
            condition_t conditions;
            for(auto const& parent : node->cpt.condition_node())
                conditions[parent] = pattern[index_node_in_graph(parent)];

            auto it = evidence.find(node);
            if(it == evidence.end())
            {
                // node is 'not' evidence
                auto const select = make_random_by_weight(probability_generator_(), node->cpt[conditions].second);
                pattern[index_node_in_graph(node)] = select;
            }
            else
            {
                // node is evidence
                w *= node->cpt[conditions].second[it->second];
                pattern[index_node_in_graph(node)] = it->second;
            }
        }

        return std::make_pair(pattern, w);
    }

    std::size_t index_node_in_graph(vertex_type const& node)
    {
        auto const& node_list = graph_.vertex_list();

        for(std::size_t i = 0; i < node_list.size(); ++i)
            if(node_list[i] == node) return i;

        throw std::runtime_error("");
    }

    // Using weight, select from random value
    // weightを使って，ランダム値からそれがどの選択値になるか判断する
    std::size_t make_random_by_weight(double const value, std::vector<double> const& weight) const
    {
        assert(0.0 <= value && value < 1.0);

        double total = 0.0;
        for(std::size_t i = 0; i < weight.size(); ++i)
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

