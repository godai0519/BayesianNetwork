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
        : filename_(), node_list_(), table_(), sampling_size_(0)
    {
    }

    sampler(std::string const& filename)
        : filename_(filename), node_list_(), table_(), sampling_size_(0)
    {
    }

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

    bool load_sample(std::vector<vertex_type> const& node_list, std::vector<element_type> const& table)
    {
        table_ = table;
        node_list_ = node_list;

        sampling_size_ = 0;
        for(auto const& p : table) sampling_size_ += p.num;

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

        // ライン数を数える
        std::size_t line_num = 0;
        {
            std::string line;
            while(std::getline(ifs, line)) ++line_num;
        }
        ifs.clear();
        ifs.seekg(0, std::ios::beg);

        // 読み込み先
        std::size_t sampling_size = 0;
        std::vector<element_type> table;
        table.reserve(line_num);

        // from file to table
        std::string line;
        while(std::getline(ifs, line))
        {
            // 空白でパース
            std::vector<std::string> select_strs;
            boost::algorithm::split(select_strs, line, boost::is_space(), boost::algorithm::token_compress_on);

            std::vector<std::size_t> select(select_strs.size());
            for(std::size_t i = 1; i < select_strs.size(); ++i)
                select[i - 1] = std::stoull(select_strs[i]);

            std::size_t const num = std::stoull(select_strs[0]);
            table.emplace_back(std::move(select), num);
            sampling_size += num;
        }

        node_list_ = node_list;
        table_.swap(table);
        sampling_size_ = sampling_size;

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

        // 指定ノード，それぞれに対してCPTを構成する
        for(auto const& node : target_nodes)
        {
            auto const parents = graph.in_vertexes(node);
            node->cpt.assign(parents, node);

            auto joint_node = parents;
            joint_node.push_back(node);

            std::unordered_map<condition_t, std::vector<std::size_t>> counter;
            for(auto const& sample : jointed_table(joint_node))
            {
                // 条件を畳み込む
                condition_t conditional;
                for(std::size_t i = 0; i < parents.size(); ++i)
                    conditional[parents[i]] = sample.select[i];

                // counterのカウントアップ
                auto it = counter.lower_bound(conditional);
                if(it == counter.end() || it->first != conditional)
                {
                    // counterに要素を追加
                    it = counter.emplace_hint(
                        it, std::piecewise_construct,
                        std::forward_as_tuple(conditional),
                        std::forward_as_tuple(node->selectable_num, 0));
                }

                // 足しあわせ
                it->second[sample.select.back()] += sample.num;
            }

            for(auto const& conditional : node->cpt.pattern())
            {
                // カウンタからデータを引きずりだし，確率1にすべき母数を算出
                auto const it = counter.find(conditional);
                auto const select_counter =
                    (it != counter.end()) ? it->second
                                          : std::vector<std::size_t>(node->selectable_num, 0);
                auto const parameter = static_cast<double>(std::accumulate(select_counter.begin(), select_counter.end(), static_cast<std::size_t>(0)));

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
        return make_cpt(graph, {target_node});
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
    std::vector<element_type> table() const
    {
        return table_;
    }

    std::vector<vertex_type> node_list() const
    {
        return node_list_;
    }

    inline std::size_t index_node(vertex_type const& node) const
    {
        for(std::size_t i = 0; i < node_list_.size(); ++i)
            if(node_list_[i] == node) return i;

        throw std::runtime_error("No sampled node");
    }

    // 周辺化テーブル
    std::vector<element_type> jointed_table(std::vector<vertex_type> const& nodes) const
    {
        // 移動先リストを作る
        std::vector<std::size_t> source_dictionary(nodes.size(), 0);
        for(std::size_t i = 0; i < nodes.size(); ++i)
        {
            for(std::size_t j = 0; j < node_list_.size(); ++j)
            {
                if(nodes[i] == node_list_[j])
                {
                    source_dictionary[i] = j;
                    break;
                }
            }
        }

        // サンプルを回して周辺化
        std::vector<element_type> joint;
        joint.reserve(table_.size());
        for(auto const& elem : table_)
        {
            // 必要なものだけ抽出
            std::vector<std::size_t> select(nodes.size());
            for(std::size_t i = 0; i < nodes.size(); ++i)
                select[i] = elem.select[source_dictionary[i]];

            // テーブルに追加されてれば加算，なければ追加
            auto const it = std::find_if(
                joint.begin(), joint.end(),
                [&select](element_type const& elem){ return elem.select == select; }
                );

            if(it == joint.end())
                joint.emplace_back(std::move(select), elem.num);
            else
                it->num += elem.num;
        }

        joint.shrink_to_fit();
        return joint;
    }

    // サンプル数に関するgetter
    std::size_t sampling_size() const
    {
        return sampling_size_;
    }

private:
    std::string filename_;
    std::vector<vertex_type> node_list_;
    std::vector<element_type> table_;
    std::size_t sampling_size_;
};

} // namespace bn

#endif // #ifndef BNI_SAMPLER_HPP
