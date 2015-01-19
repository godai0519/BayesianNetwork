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

