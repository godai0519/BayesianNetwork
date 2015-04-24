#ifndef BNI_SAMPLER_HPP
#define BNI_SAMPLER_HPP

#include <fstream>
#include <sstream>
#include <functional>
#include <vector>
#include <unordered_map>
#include <boost/optional.hpp>
#include <bayesian/graph.hpp>

namespace bn {
    
// 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
// TODO: class bpの中と重複しているので，後で整理する
void all_combination_pattern(
    std::vector<vertex_type> const& combination,
    std::function<void(condition_t const&)> const& function
    )
{
    typedef std::vector<vertex_type>::const_iterator iterator_type;
    std::function<void(iterator_type const, iterator_type const&)> recursive;

    condition_t condition;
    recursive = [&](iterator_type const it, iterator_type const& end)
    {
        if(it == end)
        {
            function(condition);
        }
        else
        {
            for(int i = 0; i < (*it)->selectable_num; ++i)
            {
                condition[*it] = i;
                recursive(it + 1, end);
            }
        }
    };

    recursive(combination.cbegin(), combination.cend());
}

class sampler {
public:
    typedef std::function<void(condition_t const&)> Handler;

    // 引数指定されたファイルからCSVデータを読み込み，1行ずつhandlerを実行
    // 引数: graph: 読み込みしたCPTを保存するグラフ構造
    //       filename: タブ区切りのCSV(要確認)
    //       node_list: ファイルに記述されている順番を指定するvector(省略可)
    //       handler: サンプルごとに実行される関数
    // 戻り値: 成功した場合のみ true
    bool load_sample(graph_t const& graph, std::string const& filename, Handler handler) const;
    bool load_sample(graph_t const& graph, std::string const& filename, std::vector<vertex_type> const& node_list, Handler handler) const;
    
    // 引数指定されたファイルからCSVデータを読み込み，CSVデータに従ってCPTを作成する
    // 引数: graph: 読み込みしたCPTを保存するグラフ構造
    //       filename: タブ区切りのCSV(要確認)
    //       node_list: ファイルに記述されている順番を指定するvector(省略可)
    // 戻り値: 成功した場合のみ true
    bool load_cpt(graph_t const& graph, std::string const& filename) const;
    bool load_cpt(graph_t const& graph, std::string const& filename, std::vector<vertex_type> const& node_list) const;

    // Sampling数のgetter
    // sampling数は，頂点/辺の追加削除によって信用できなくなることに注意すること
    boost::optional<std::size_t> const& sampling_size() const;

private:
    mutable boost::optional<std::size_t> sampling_size_;
};


bool sampler::load_sample(graph_t const& graph, std::string const& filename, Handler handler) const
{
    return load_sample(graph, filename, graph.vertex_list(), handler);
}

bool sampler::load_sample(graph_t const& graph, std::string const& filename, std::vector<vertex_type> const& node_list, Handler handler) const
{
    // ファイルオープン
    std::ifstream ifs(filename);
    if(!ifs.is_open()) return false;

    // 1行ずつサンプリング
    std::string line_str;
    std::size_t sample_num = 0;
    while(std::getline(ifs, line_str))
    {
        condition_t sample;
        std::transform(
            node_list.cbegin(), node_list.cend(),
            std::istream_iterator<std::string>(std::istringstream(line_str) >> std::skipws),
            std::inserter(sample, sample.begin()),
            [](vertex_type const& node, std::string const& str){ return std::make_pair(node, std::stoi(str)); }
            );

        handler(sample);
        ++sample_num;
    }

    sampling_size_ = sample_num;
    return true;
}

bool sampler::load_cpt(graph_t const& graph, std::string const& filename) const
{
    return load_cpt(graph, filename, graph.vertex_list());
}

bool sampler::load_cpt(graph_t const& graph, std::string const& filename, std::vector<vertex_type> const& node_list) const
{
     // CPTの初期化
    for(auto const& node : graph.vertex_list())
    {
        auto const in_vertex = graph.in_vertexes(node);
        node->cpt.assign(in_vertex, node);
    }

    // カウンタの初期化
    std::unordered_map<vertex_type, std::unordered_map<condition_t, std::vector<int>>> match_counter;
    std::unordered_map<vertex_type, std::unordered_map<condition_t, int>> counter;
    for(auto const& node : graph.vertex_list())
    {
        auto const in_vertex = graph.in_vertexes(node);
        all_combination_pattern(
            in_vertex,
            [&node, &counter, &match_counter](condition_t const& condition)
            {
                counter[node][condition] = 0;
                match_counter[node][condition].assign(node->selectable_num, 0);
            });
    }
    
    // 1行ずつサンプリング
    load_sample(
        graph, filename,
        [&graph, &counter, &match_counter](condition_t const& sample)
        {
            for(auto const& node : graph.vertex_list())
            {
                auto const in_vertex = graph.in_vertexes(node);
                all_combination_pattern(
                    in_vertex,
                    [&sample, &node, &counter, &match_counter](condition_t const& condition) -> void
                    {
                        // 条件に合うかどうか
                        bool const meet = std::all_of(
                            condition.cbegin(), condition.cend(),
                            [&sample](condition_t::value_type const& p) -> bool{ return sample.at(p.first) == p.second; });

                        if(meet)
                        {
                            ++counter[node][condition];
                            ++match_counter[node][condition][sample.at(node)];
                        }
                    });
            }
            return;
        });

    // サンプリングの正規化
    for(auto const& node : graph.vertex_list())
    {
        auto const in_vertex = graph.in_vertexes(node);
        all_combination_pattern(
            in_vertex,
            [&node, &counter, &match_counter](condition_t const& condition)
            {
                for(std::size_t i = 0; i < node->selectable_num; ++i)
                {
                    if(counter[node][condition] != 0)
                        node->cpt[condition].second[i] = (double)match_counter[node][condition][i] / counter[node][condition];
                }
            });
    }

    return true;
}


boost::optional<std::size_t> const& sampler::sampling_size() const
{
    return this->sampling_size_;
}


} // namespace bn

#endif // #ifndef BNI_SAMPLER_HPP
