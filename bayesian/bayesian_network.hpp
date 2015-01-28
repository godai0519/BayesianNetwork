#ifndef BNI_BAYESIAN_NETWORK_HPP
#define BNI_BAYESIAN_NETWORK_HPP

#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include "graph.hpp"

namespace bn {

template<class NodeType>
class bayesian_network {
public:
    // グラフをセットしないctor
    explicit bayesian_network();

    // 引数指定されたファイルからデータを読み込み，データ系列をメンバ変数data_に格納する
    // 引数: filename: タブ区切りのCSV
    //       node_list: ファイルに記述されている順番を指定するvector
    // 成功した場合，trueが返される
    bool load_data(std::string const& filename, std::vector<vertex_type> const& node_list);

    // 事前に読みだしたCSVデータに従ってCPTを作成する
    // 引数: 読み込みしたCPTを保存するgraph_t
    // 成功した場合，trueが返される
    bool load_cpt(graph_t const& graph);

    // load_cptされたrawデータ系列を返す
    std::vector<condition_t> data() const;

    // 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
    // TODO: class bpの中と重複しているので，後で整理する
    void all_combination_pattern(
        std::vector<vertex_type> const& combination,
        std::function<void(condition_t const&)> const& function
        );

private:
    std::vector<condition_t> data_;
};

template<class NodeType>
bayesian_network<NodeType>::bayesian_network()
{
}

template<class NodeType>
bool bayesian_network<NodeType>::load_data(std::string const& filename, std::vector<vertex_type> const& node_list)
{
    // オープンして確認
    std::ifstream ifs(filename);
    if(!ifs.is_open()) return false;

    // 読み込み
    data_.clear();
    std::string line;
    while(std::getline(ifs, line))
    {
        condition_t cond;
        std::istringstream iss(line);
        std::transform(
            node_list.cbegin(), node_list.cend(),
            std::istream_iterator<std::string>(iss),
            std::inserter(cond, cond.begin()),
            [](vertex_type const& node, std::string const& str){ return std::make_pair(node, std::stoi(str)); }
            );
        data_.push_back(std::move(cond));
    }

    return true;
}

template<class NodeType>
bool bayesian_network<NodeType>::load_cpt(graph_t const& graph)
{
    // データが読み込みされていない可能性がある
    if(data_.size() == 0) return false;

    // 読み込みを行う
    for(auto const& node : graph_.vertex_list())
    {
        auto const in_vertex = graph_.in_vertexs(node);
        node->cpt.assign(in_vertex, node); // CPTの初期化

        // 親ノードのすべての組み合わせに対し，ループを回し，該当するサンプルがあった場合は数え上げる
        all_combination_pattern(
            in_vertex,
            [this, &node](condition_t const& condition)
            {
                int count = 0;
                auto& cpt = node->cpt[condition].second;
                for(auto const& sample : data_)
                {
                    // 条件に合うかどうか
                    bool const meet = std::all_of(
                        condition.cbegin(), condition.cend(),
                        [&sample](condition_t::value_type const& p){ return sample.at(p.first) == p.second; });

                    if(meet){
                        ++count;
                        cpt[sample.at(node)] += 1.0;
                    }
                }

                for(std::size_t i = 0; i < cpt.size(); ++i)
                {
                    cpt[i] /= count;
                }
            });
    }

    return true;
}

template<class NodeType>
std::vector<condition_t> bayesian_network<NodeType>::data() const
{
    return data_;
}

template<class NodeType>
void bayesian_network<NodeType>::all_combination_pattern(
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

} // namespace bn

#endif // #ifndef BNI_BAYESIAN_NETWORK_HPP

