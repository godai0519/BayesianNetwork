#ifndef BNI_EVALUATION_TRANSINFORMATION_HPP
#define BNI_EVALUATION_TRANSINFORMATION_HPP

#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>

namespace bn {
namespace evaluation {

struct entropy {
    // samplingから得られるテーブルを用いて，variablesの同時entropyを計算する
    // sampling: load_sample済のsampler
    // variables: 調べるノードのリスト
    double operator() (sampler const& sampling, std::vector<vertex_type> const& variables) const
    {
        // 元データから厳密な確率を出す(条件付き確率ではない)
        // CPTを親ノードを行列乗算して算出することも可能だが，どうだろうか
        std::unordered_map<condition_t, std::size_t> table;
        condition_t cond;
        for(auto const& sample : sampling.table())
        {
            for(auto const& variable : variables)
                cond[variable] = sample.first.at(variable);

            // 有れば加算，なければ作成
            auto it = table.lower_bound(cond);
            if(it != table.end() && !(table.key_comp()(cond, it->first)))
                it->second += sample.second;
            else
                table.insert(it, std::make_pair(cond, sample.second));
        }

        // すべての確率に底2のエントロピー計算を行う
        double entropy = 0.0;
        for(auto const& pattern : table)
        {
            auto const probability = static_cast<double>(pattern.second) / sampling.sampling_size();
            entropy -= probability * std::log2(probability);
        }
        return entropy;
    }

    // samplingから得られるテーブルを用いて，variablesの同時entropyを計算する
    // sampling: load_sample済のsampler
    // variable: 調べるノード
    double operator() (sampler const& sampling, vertex_type const& variable) const
    {
        std::vector<vertex_type> variables = {variable};
        return (*this)(sampling, variables);
    }
};

struct mutual_information {
    // samplingから得られるテーブルを用いて，xとyの間の相互情報量を計算する
    // sampling: load_sample済のsampler
    // x, y: 調べるノード
    double operator() (sampler const& sampling, vertex_type const& x, vertex_type const& y) const
    {
        // Initialize table
        std::vector<double> x_counter(x->selectable_num, 0.0);
        std::vector<double> y_counter(y->selectable_num, 0.0);
        std::vector<std::vector<double>> joint_counter(
            x->selectable_num, std::vector<double>(y->selectable_num, 0.0)
            );

        // Count up
        for(auto const& sample : sampling.table())
        {
            auto const x_value = sample.first.at(x);
            auto const y_value = sample.first.at(y);
            auto const probability = static_cast<double>(sample.second) / sampling.sampling_size();
            x_counter[x_value] += probability;
            y_counter[y_value] += probability;
            joint_counter[x_value][y_value] += probability;
        }

        // Calculate
        double mi_value = 0.0;
        for(std::size_t x_value = 0; x_value < x->selectable_num; ++x_value)
        {
            for(std::size_t y_value = 0; y_value < y->selectable_num; ++y_value)
            {
                auto const multiply_probability = x_counter[x_value] * y_counter[y_value];
                auto const joint_probability = joint_counter[x_value][y_value];

                if(joint_probability == 0) break; // 0 * log0 -> 0
                else mi_value += joint_probability * std::log2(joint_probability / multiply_probability);
            }
        }

        return mi_value;
    }
};

} // namespace evaluation
} // namespace bn

#endif // #ifndef BNI_EVALUATION_TRANSINFORMATION_HPP
