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

    // 引数のグラフをmove or copyで格納するctor
    template<class T>
    explicit bayesian_network(T && graph);

    // グラフがセットされているかどうか
    bool is_graph() const;
    
    // 引数のグラフをmove or copyで格納する
    template<class T>
    void set_graph(T && graph);
   
    // 登録されたグラフを削除する
    void reset_graph();

    // 現在格納されているグラフを返す
    // is_graph() == falseのときは，空のグラフを返す
    graph_t graph() const;

    // 引数指定されたファイルから，データを読み込みCPTを生成する
    // 成功した場合，trueが返される
    // CPTはグラフの中，rawデータ系列はメンバ変数data_に格納される
    bool load_cpt(std::string const& filename, std::vector<vertex_type> const& node_list);

    // load_cptされたrawデータ系列を返す
    std::vector<condition_t> data() const;

    // 与えられた確率変数全ての組み合わせに対し，functionを実行するというインターフェースを提供する
    // TODO: class bpの中と重複しているので，後で整理する
    void all_combination_pattern(
        std::vector<vertex_type> const& combination,
        std::function<void(condition_t const&)> const& function
        );

/*
    template<class Func>
    auto apply(Func const& f) -> decltype(f(graph_));
*/
private:
    bool is_graph_ = false;
    graph_t graph_;
    std::vector<condition_t> data_;
};

template<class NodeType>
bayesian_network<NodeType>::bayesian_network()
{
}

template<class NodeType> template<class T>
bayesian_network<NodeType>::bayesian_network(T && graph)
{
    set_graph(std::forward<T>(graph));
}

template<class NodeType>
bool bayesian_network<NodeType>::is_graph() const
{
    return is_graph_;
}

template<class NodeType> template<class T>
void bayesian_network<NodeType>::set_graph(T && graph)
{
    is_graph_ = true;
    graph_ = std::forward<T>(graph);
}

template<class NodeType>
void bayesian_network<NodeType>::reset_graph()
{
    is_graph_ = false;
    graph_.swap(graph_t());
}

template<class NodeType>
graph_t bayesian_network<NodeType>::graph() const
{
    if(is_graph()) return graph_;
    else return graph_t();
}

template<class NodeType>
bool bayesian_network<NodeType>::load_cpt(std::string const& filename, std::vector<vertex_type> const& node_list)
{
    // グラフがセットされていなかった場合はエラー
    if(!is_graph()) return false;

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
    
    for(auto const& node : node_list)
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

/*
template<class NodeType> template<class Func>
auto bayesian_network<NodeType>::apply(Func const& f) -> decltype(f(graph_))
{
    return f(graph);
}
*/

} // namespace bn

#endif // #ifndef BNI_BAYESIAN_NETWORK_HPP

