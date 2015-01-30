#ifndef BNI_LIKELIHOOD_WEIGHTING_HPP
#define BNI_LIKELIHOOD_WEIGHTING_HPP

#include <random>
#include "graph.hpp"

namespace bn {

class likelihood_weighting {
public:
    typedef std::unordered_map<vertex_type, int> evidence_list;
    typedef std::unordered_map<vertex_type, int> pattern_list;
    typedef std::unordered_map<vertex_type, matrix_type> return_type;

    explicit likelihood_weighting(graph_t const& graph);
    virtual ~likelihood_weighting() = default;

    // Run: Likelihood Weighting
    return_type operator() (evidence_list const& evidence, std::uint64_t const sample_num = 10000);

private:
    // Make a sample according to evidence node
    // evidenceノードに従って，1サンプルを作成する
    std::pair<pattern_list, double> weighted_sample(evidence_list const& evidence);
    
    // Using weight, select from random value
    // weightを使って，ランダム値からそれがどの選択値になるか判断する
    int make_random_by_weight(double const value, std::vector<double> const& weight) const;

    // Normalization Σ[i,j](a_ij) = 1
    // (重複)
    matrix_type& normalize(matrix_type& target) const;
    
    // 乱数生成器
    class probability_generator {
    public:
        probability_generator();
        double operator() ();
    private:
        std::unique_ptr<std::mt19937> engine_;
        std::uniform_real_distribution<double> distribution_;
    };

    graph_t const graph_;
    probability_generator probability_generator_;
};

} // namespace bn

#endif // #ifndef BNI_LIKELIHOOD_WEIGHTING_HPP

