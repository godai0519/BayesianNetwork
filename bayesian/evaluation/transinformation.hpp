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
        for(auto const& sample : sampling.table())
        {
            condition_t cond;
            for(auto const& variable : variables)
                cond[variable] = sample.first.at(variable);

            // 有れば加算，なければ作成
            auto it = table.find(cond);
            if(it == table.end()) table.emplace(cond, sample.second);
            else                  it->second += sample.second;
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
        entropy ent;
        return ent(sampling, x) + ent(sampling, y) - ent(sampling, {x, y});
    }

    // samplingから得られるテーブルを用いて，xとyの間の相互情報量を計算する(xとy単体のentropyが計算済の場合)
    // sampling: load_sample済のsampler
    // x, y: 調べるノード
    // x_ent, y_ent: 調べるノードのエントロピー
    template<class T>
    double operator() (sampler const& sampling, vertex_type const& x, T const x_ent, vertex_type const& y, T const y_ent) const
    {
        entropy ent;
        return x_ent + y_ent - ent(sampling, {x, y});
    }

    // samplingから得られるテーブルを用いて，xとyの間の相互情報量を計算する(xとyとxyのentropyが計算済の場合)
    // x_ent, y_ent: 調べるノードのエントロピー
    // xy_ent: 調べるノード2つの同時エントロピー
    template<class T>
    double operator() (T const x_ent, T const y_ent, T const xy_ent) const
    {
        return x_ent + y_ent - xy_ent;
    }
};

} // namespace evaluation
} // namespace bn

#endif // #ifndef BNI_EVALUATION_TRANSINFORMATION_HPP
