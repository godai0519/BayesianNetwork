#ifndef BNI_SAMPLER_HPP
#define BNI_SAMPLER_HPP

#include <fstream>
#include <sstream>
#include <functional>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <boost/optional.hpp>
#include <boost/algorithm/string.hpp>
#include <bayesian/graph.hpp>
#include <bayesian/utility.hpp>
#include <bayesian/hash.hpp>

namespace bn {

class sampler {
public:
    sampler()
        : filename_(), table_(), sampling_size_(0)
    {
    }

    sampler(std::string const& filename)
        : filename_(filename), table_(), sampling_size_(0)
    {
    }

    bool load_sample(std::unordered_map<condition_t, std::size_t> const& table)
    {
        table_ = table;

        sampling_size_ = 0;
        for(auto const& p : table) sampling_size_ += p.second;

        return true;
    }

    // 事前に設定したファイルからCSVデータを読み込み，一般化したテーブルを作る
    // 引数: node_list: ファイルに記述されている順番を指定するvector
    // 戻り値: 成功した場合のみ true
    bool load_sample(std::vector<vertex_type> const& node_list)
    {
        // ファイルオープン
        std::ifstream ifs(filename_);
        if(!ifs.is_open()) return false;

        // 読み込み先
        std::size_t sampling_size = 0;
        std::unordered_map<condition_t, std::size_t> table;

        // 1行ずつサンプリング
        std::string line_str;
        condition_t sample;
        while(std::getline(ifs, line_str))
        {
            std::vector<std::string> line;
            boost::algorithm::split(line, line_str, boost::is_space(), boost::algorithm::token_compress_on);

            auto select_it = line.cbegin() + 1;
            for(auto it = node_list.cbegin(); it != node_list.cend();)
                sample[*it++] = std::stoi(*select_it++);

            auto const sample_num = std::stoi(line[0]);

            auto const it = table.find(sample);
            if(it == table.cend()) table[sample] = sample_num;
            else (it->second) += sample_num;

            sampling_size += sample_num;
        }

        sampling_size_ = sampling_size;
        table_ = std::move(table);
        return true;
    }

	// load_sampleによって生成したテーブルから，CPTを計算し，graphに格納する
	// 引数: graph: 読み込みしたCPTを保存するグラフ構造
	// 戻り値: 成功した場合のみ true
    bool make_cpt(graph_t const& graph) const
    {
        return make_cpt(graph, graph.vertex_list());
    }

	// load_sampleによって生成したテーブルから，CPTを計算し，graphに格納する
	// 引数: graph: 読み込みしたCPTを保存するグラフ構造
    //       target_nodes: 読み込みするノード群
	// 戻り値: 成功した場合のみ true
    bool make_cpt(graph_t const& graph, std::vector<vertex_type> const& target_nodes) const
    {
        // sampleが読み込まれているかどうか
        if(sampling_size() == 0) return false;

        // カウンタ
        std::unordered_map<
            std::pair<vertex_type, condition_t>,
            std::vector<std::size_t>,
            bn::hash<std::pair<vertex_type, condition_t>>
        > counter;

        // カウント
        for(auto const& sample : table_)
        {
            for(auto const& node : target_nodes)
            {
                auto const parents = graph.in_vertexes(node);

                // 条件を畳み込む
                condition_t conditional;
                for(auto const& parent : parents) conditional[parent] = sample.first.at(parent);

                // counterのカウントアップ
                auto container_first = std::make_pair(node, std::move(conditional));
                auto it = counter.lower_bound(container_first);
                if(it == counter.end() || it->first != container_first)
                {
                    // counterに要素を追加
                    it = counter.emplace_hint(
                        it, std::piecewise_construct,
                        std::forward_as_tuple(container_first),
                        std::forward_as_tuple(node->selectable_num, 0));
                }
                
                // 足しあわせ
                it->second[sample.first.at(node)] += sample.second;
            }
        }

        // カウント値を用い，CPTを算出
        for(auto const& node : target_nodes)
        {
            // CPTの初期化
            auto const in_vertex = graph.in_vertexes(node);
            node->cpt.assign(in_vertex, node);

            for(auto const& conditional : node->cpt.pattern())
            {
                // カウンタからデータを引きずりだし，確率1にすべき母数を算出
                auto const key = std::make_pair(node, conditional);
                auto const it = counter.find(key);
                auto const select_counter =
                    (it != counter.end()) ? counter.at(std::make_pair(node, conditional))
                                          : std::vector<std::size_t>(node->selectable_num, 0);

                double const parameter = static_cast<double>(std::accumulate(select_counter.begin(), select_counter.end(), static_cast<std::size_t>(0)));

                // 本題のCPT作成
                std::vector<double> cpt(node->selectable_num);
                if(parameter == 0)
                {
                    // 同じ値を与える
                    for(std::size_t i = 0; i < node->selectable_num; ++i)
                        cpt[i] = static_cast<double>(1.0) / node->selectable_num;
                }
                else
                {
                    // CPTを正規化する
                    for(std::size_t i = 0; i < node->selectable_num; ++i)
                        cpt[i] = static_cast<double>(select_counter.at(i)) / parameter;
                }

                node->cpt[conditional].second = std::move(cpt);
            }
        }

        return true;
    }

	// load_sampleによって生成したテーブルから，CPTを計算し，graphに格納する
	// 引数: graph: 読み込みしたCPTを保存するグラフ構造
    //       target_node: 読み込みするノード
	// 戻り値: 成功した場合のみ true
    bool make_cpt(graph_t const& graph, vertex_type const& target_node) const
    {
        // sampleが読み込まれているかどうか
        if(sampling_size() == 0) return false;

        // カウンタ
        std::unordered_map<
            condition_t,
            std::vector<std::size_t>
        > counter;
        
        // CPTの初期化
        auto const parents = graph.in_vertexes(target_node);
        target_node->cpt.assign(parents, target_node);

        // カウント
        for(auto const& sample : table_)
        {

            // 条件を畳み込む
            condition_t conditional;
            for(auto const& parent : parents) conditional[parent] = sample.first.at(parent);

            // counterのカウントアップ
            auto it = counter.lower_bound(conditional);
            if(it == counter.end() || it->first != conditional)
            {
                // counterに要素を追加
                it = counter.emplace_hint(
                    it, std::piecewise_construct,
                    std::forward_as_tuple(conditional),
                    std::forward_as_tuple(target_node->selectable_num, 0));
            }
                
            // 足しあわせ
            it->second[sample.first.at(target_node)] += sample.second;
        }

        for(auto const& conditional : target_node->cpt.pattern())
        {
            // カウンタからデータを引きずりだし，確率1にすべき母数を算出
            auto const it = counter.find(conditional);
            auto const select_counter =
                (it != counter.end()) ? it->second
                                      : std::vector<std::size_t>(target_node->selectable_num, 0);
            auto const parameter = static_cast<double>(std::accumulate(select_counter.begin(), select_counter.end(), static_cast<std::size_t>(0)));

            // 本題のCPT作成
            std::vector<double> cpt(target_node->selectable_num);
            if(parameter == 0)
            {
                // 同じ値を与える
                for(std::size_t i = 0; i < target_node->selectable_num; ++i)
                    cpt[i] = static_cast<double>(1.0) / target_node->selectable_num;
            }
            else
            {
                // CPTを正規化する
                for(std::size_t i = 0; i < target_node->selectable_num; ++i)
                    cpt[i] = static_cast<double>(select_counter.at(i)) / parameter;
            }

            target_node->cpt[conditional].second = std::move(cpt);
        }

        return true;
    }

    // 読み込み対象のファイルを指定するgetter
    std::string filename() const
    {
        return filename_;
    }

    // 読み込み対象のファイルを指定するsetter
    void set_filename(std::string const& filename)
    {
        filename_ = filename;
        sampling_size_ = 0;
        table_.clear();
    }

	// テーブルに関するgetter
	std::unordered_map<condition_t, std::size_t> table() const
    {
        return table_;
    }

	// サンプル数に関するgetter
    std::size_t sampling_size() const
    {
        return sampling_size_;
    }

private:
    std::string filename_;
    std::unordered_map<condition_t, std::size_t> table_;
    std::size_t sampling_size_;
};

} // namespace bn

#endif // #ifndef BNI_SAMPLER_HPP
