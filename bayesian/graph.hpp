#ifndef BNI_GRAPH_HPP
#define BNI_GRAPH_HPP

#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>
#include <cstdint>
#include <memory>
#include <unordered_map>

namespace bn {

struct vertex_t;
struct edge_t;
typedef std::shared_ptr<vertex_t> vertex_type;
typedef std::shared_ptr<edge_t>   edge_type;
typedef std::unordered_map<vertex_type, int> condition_t;

} // namespace bn

namespace std {

// from Boost
template<class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// ぐぬぬ．規格的には特殊化はstd実装許可されてる(晶さん)
#ifdef _MSC_VER
template<> struct hash<bn::condition_t> : public unary_function<bn::condition_t, size_t> {
    inline std::size_t operator()(bn::condition_t const& cond) const
#else
template<> struct hash<bn::condition_t> : public __hash_base<std::size_t, bn::condition_t> {
    inline std::size_t operator()(bn::condition_t const& cond) const noexcept
#endif
    {
        std::size_t value = 0;
        for(auto const& one : cond)
        {
            value ^= std::hash<bn::vertex_type>()(one.first);
            value ^= std::hash<int>()(one.second);
        }
        return value;
    }
};

} // namespace std

namespace bn {

// Conditional Probability Table (CPT) 表現クラス
// 主に vertex_t が保持する
class cpt_t {
public:
    typedef std::unordered_map<condition_t, std::vector<double>> table_type;

    explicit cpt_t() = default;
    explicit cpt_t(std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node)
    {
        assign(parent_nodes, target_node);
    }

    // Select Parent Node (Conditional) and Target Node
    // Use them to make CPT
    inline void assign(std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node);

    // 一部の条件を元にそれに該当するデータを抽出してリストとして返す
    // 空条件を渡すと，実質的にはtable()と同じ
    table_type filter(condition_t const& cond)
    {
        table_type filtered;
        std::copy_if(
            table_.begin(), table_.end(), std::inserter(filtered, filtered.begin()),
            [&cond](std::pair<condition_t, std::vector<double>> const& target) -> bool
            {
                for(auto it = cond.begin(); it != cond.end(); ++it)
                {
                    auto result = target.first.find(it->first);
                    if(result == target.first.end() || result->second != it->second)
                    {
                        return false;
                    }
                }

                return true;
            });

        return filtered;
    }

    // 条件として使用できるノードリストを返す
    std::vector<vertex_type> condition_node()
    {
        return parents_;
    }

    // operator[]でアクセス可能な全てのパターンを返却する
    std::vector<condition_t> pattern()
    {
        std::vector<condition_t> pat;
        pat.reserve(table_.size());

        for(auto const& elem : table_)
            pat.push_back(elem.first);

        return pat;
    }

    // 条件が完全一致した確率vectorを返す
    // std::pairのfirstが検索成功したかをのせる
    // firstがtrueのとき，secondには実体への参照が格納される
    // firstがfalseのときのsecondについては未定義
    std::pair<bool, std::vector<double>&> operator[] (condition_t const cond)
    {
        auto result = table_.find(cond);
        if(result == table_.end())
        {
            std::vector<double> dummy;
            return std::pair<bool, std::vector<double>&>(false, std::ref(dummy));
        }
        else
        {
            return std::pair<bool, std::vector<double>&>(true, std::ref(result->second));
        }
    }

    // 条件が完全一致した確率vectorを返す (const版)
    // std::pairのfirstが検索成功したかをのせる
    // firstがtrueのとき，secondには実体への参照が格納される
    // firstがfalseのときのsecondについては未定義
    std::pair<bool, std::vector<double> const&> operator[] (condition_t const cond) const
    {
        auto result = table_.find(cond);
        if(result == table_.end())
        {
            std::vector<double> dummy;
            return std::pair<bool, std::vector<double> const&>(false, std::ref(dummy));
        }
        else
        {
            return std::pair<bool, std::vector<double> const&>(true, std::ref(result->second));
        }
    }

private:
    inline void assign_impl(table_type& new_table, condition_t cond, std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node, std::size_t const n) const;

    std::vector<vertex_type> parents_;
    table_type table_;
};

// 頂点を示す(note_tに変更することを検討)
struct vertex_t {
    int id;
    std::size_t selectable_num = 0; // 取りうる値の数
    cpt_t cpt;
};

struct edge_t {}; // tag

struct database_t {
    std::string graph_name;
    std::unordered_map<std::size_t, std::string> node_name;
    std::unordered_map<std::size_t, std::vector<std::string>> options_name;
};

// ネットワーク構造を表現する
// 頂点・辺オブジェクトを入力とし，対応するオブジェクトを返却する
class graph_t {
public:
    graph_t() = default;
    virtual ~graph_t() = default;

    // copy/move ctor
    graph_t(graph_t const& other) = default;
    graph_t(graph_t && other)
    {
        other.swap(*this);
    }

    graph_t& operator=(graph_t const& rhs)
    {
        graph_t(rhs).swap(*this);
        return *this;
    }
    graph_t& operator=(graph_t&& rhs)
    {
        rhs.swap(*this);
        return *this;
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    inline void swap(graph_t& other)
#else
    void inline swap(graph_t& other) noexcept
#endif
    {
        std::swap(vertex_list_, other.vertex_list_);
        std::swap(edge_list_, other.edge_list_);
        std::swap(adjacent_list_, other.adjacent_list_);
    }


#if defined(_MSC_VER) && _MSC_VER < 1900
    friend inline void swap(graph_t& lhs, graph_t& rhs);
#else
    friend inline void swap(graph_t& lhs, graph_t& rhs) noexcept;
#endif

    std::vector<vertex_type> const& vertex_list() const
    {
        return vertex_list_;
    }

    std::vector<edge_type> const& edge_list() const
    {
        return edge_list_;
    }

    graph_t clone() const
    {
        graph_t cloned_graph;

        auto const& vertex_list = this->vertex_list();
        auto const& cloned_vertex_list = cloned_graph.vertex_list();

        // 頂点のコピー
        for(auto const& vertex : vertex_list)
        {
            auto cloned_vertex = cloned_graph.add_vertex();
            *cloned_vertex = *vertex;
        }

        // 辺のコピー
        for(auto const& edge : this->edge_list())
        {
            auto const source = this->index_search(this->source(edge));
            auto const target = this->index_search(this->target(edge));
            cloned_graph.add_edge(cloned_vertex_list[source], cloned_vertex_list[target]);
        }

        return cloned_graph;
    }

    // 頂点を生成し，そのshared_ptrを返す
    // 必ず成功する
    vertex_type add_vertex()
    {
        // 登録
        auto v = std::make_shared<vertex_t>();
        vertex_list_.push_back(v);

        // adjacent_list_を1行1段増やす
        auto const new_size = vertex_list_.size();
        for(auto& line : adjacent_list_) line.resize(new_size, nullptr); // 既存のリサイズ
        adjacent_list_.emplace_back(new_size, nullptr); // 新規追加

        return v;
    }

    // 引数を元に辺を生成し，そのshared_ptrを返す
    // 既に辺が存在していた場合やDAGが成立しなくなる場合は，nullptrなshared_ptrを即座に返却する
    // 頂点が存在していなかった場合も同様
    edge_type add_edge(vertex_type const& from, vertex_type const& to)
    {
        if(is_able_trace(to, from))
        {
            // これから張ろうとしている辺の終点から，始点に戻る経路があるとき失敗させる
            // (閉路になるため，DAGを満たさなくなるから)
            return nullptr;
        }

        auto const index_from = index_search(from);
        auto const index_to   = index_search(to);

        // indexの溢れがあるかどうか
        if(index_from >= vertex_list_.size() || index_to >= vertex_list_.size()) return nullptr;

        // 辺の存在を確認
        if(adjacent_list_[index_from][index_to]) return nullptr;

        // 登録
        auto e = std::make_shared<edge_t>();
        edge_list_.push_back(e);
        adjacent_list_[index_from][index_to] = e;
        return e;
    }

    // 引数を元に一致する頂点を削除する
    // その際，その頂点に入出辺も同時に削除する
    // 削除に成功した場合はtrueを返却する
    bool erase_vertex(vertex_type const& v)
    {
        auto const index = index_search(v);
        if(index >= vertex_list_.size()) return false;

        vertex_list_.erase(vertex_list_.begin() + index);
        adjacent_list_.erase(adjacent_list_.begin() + index);
        for(auto& line : adjacent_list_) line.erase(line.begin() + index);
        return true;
    }

    // 引数を元に一致する辺を削除する
    // 削除に成功した場合はtrueを返却する
    bool erase_edge(edge_type const& e)
    {
        auto index = edge_search(e);
        if(index.first >= vertex_list_.size() || index.second >= vertex_list_.size()) return false;

        edge_list_.erase(std::remove(edge_list_.begin(), edge_list_.end(), e), edge_list_.end());
        adjacent_list_[index.first][index.second] = nullptr;
        return true;
    }

    // 全ての頂点を削除する(付随的にerase_all_edgeも起こる)
    bool erase_all_vertex()
    {
        vertex_list_.clear();
        edge_list_.clear();
        adjacent_list_.clear();
        return true;
    }

    // 全ての辺を削除する
    bool erase_all_edge()
    {
        edge_list_.clear();
        adjacent_list_.assign(vertex_list_.size(), std::vector<edge_type>(vertex_list_.size(), nullptr));
        return true;
    }

    // 引数で指定された辺を逆向きにする
    // 渡された辺は，成功・失敗に限らずに無効化させられる
    // 失敗した場合はnullptr，成功した場合は変更された辺を返す
    edge_type change_edge_direction(edge_type const& e)
    {
        auto const to   = target(e);
        auto const from = source(e);

        // 削除
        if(!erase_edge(e))
            return nullptr;

        // 追加
        if(auto const new_edge = add_edge(to, from))
        {
            return new_edge;
        }
        else
        {
            //ダメだったら戻す
            add_edge(from, to);
            return nullptr;
        }
    }

    // 引数の頂点から出て行く辺を列挙する
    std::vector<edge_type> out_edges(vertex_type const& from) const
    {
        auto const index = index_search(from);

        std::vector<edge_type> ret;
        for(std::size_t i = 0; i < vertex_list_.size(); ++i)
        {
            if(adjacent_list_[index][i])
            {
                ret.push_back(adjacent_list_[index][i]);
            }
        }
        return ret;
    }

    // 引数の頂点から出て行く隣接頂点を列挙する
    std::vector<vertex_type> out_vertexes(vertex_type const& from) const
    {
        std::vector<vertex_type> ret;
        for(auto const& edge : out_edges(from))
        {
            ret.push_back(target(edge));
        }
        return ret;
    }

    // 引数の頂点へ入っていく辺を列挙する
    std::vector<edge_type> in_edges(vertex_type const& to) const
    {
        auto const index = index_search(to);

        std::vector<edge_type> ret;
        for(std::size_t i = 0; i < vertex_list_.size(); ++i)
        {
            if(adjacent_list_[i][index])
            {
                ret.push_back(adjacent_list_[i][index]);
            }
        }
        return ret;
    }

    // 引数の頂点へ入っていく隣接頂点を列挙する
    std::vector<vertex_type> in_vertexes(vertex_type const& to) const
    {
        std::vector<vertex_type> ret;
        for(auto const& edge : in_edges(to))
        {
            ret.push_back(source(edge));
        }
        return ret;
    }

    // 辺の先を探す
    vertex_type source(edge_type const& edge) const
    {
        auto const index = edge_search(edge);
        auto const tmp = std::numeric_limits<std::size_t>::max();
        if(index.first == tmp || index.second == tmp) return nullptr;

        return vertex_list_[index.first];
    }

    // 辺の元を探す
    vertex_type target(edge_type const& edge) const
    {
        auto const index = edge_search(edge);
        auto const tmp = std::numeric_limits<std::size_t>::max();
        if(index.first == tmp || index.second == tmp) return nullptr;

        return vertex_list_[index.second];
    }

    // fromの頂点からtoの頂点が辿れるならばtrue，それ以外ならfalseを返す
    // graphはDAGであることを前提とする
    bool is_able_trace(vertex_type const& from, vertex_type const& to) const
    {
        if(from == to) return true;

        auto const children = out_vertexes(from);
        auto const result = std::any_of(
            children.cbegin(), children.cend(),
            [this, &to](vertex_type const& next)
            {
                return is_able_trace(next, to);
            });
        return result;
    }

    //トポロジカルソートを使用して，DAG判定
    // http://ja.wikipedia.org/wiki/トポロジカルソート
    //bool is_dag();

private:
    // vertex_typeから，vertex_list_内のindexを見つける
    std::size_t index_search(vertex_type const& vertex) const
    {
        auto const iterator = std::find(vertex_list_.cbegin(), vertex_list_.cend(), vertex);
        return std::distance(vertex_list_.cbegin(), iterator);
    }

    // edge_typeからadjacent_list_内の座標を見つける
    // 該当しないときは
    // {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()}
    std::pair<std::size_t, std::size_t> edge_search(edge_type const& edge) const
    {
        for(std::size_t i = 0; i < vertex_list_.size(); ++i)
        {
            for(std::size_t j = 0; j < vertex_list_.size(); ++j)
            {
                if(adjacent_list_.at(i).at(j) == edge)
                {
                    return std::make_pair(i, j);
                }
            }
        }

        auto const tmp = std::numeric_limits<std::size_t>::max();
        return std::make_pair(tmp, tmp);
    }

    std::vector<vertex_type> vertex_list_;
    std::vector<edge_type>   edge_list_;
    std::vector<std::vector<edge_type>> adjacent_list_;
};

// inline Trick
// cpt_tをtemplate化すれば不要
inline void cpt_t::assign(std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node)
{
    std::size_t table_size = 1;
    for(auto const& node : parent_nodes)
    {
        table_size *= node->selectable_num;
    }

    condition_t cond;
    table_type new_table;
    new_table.reserve(table_size);
    assign_impl(new_table, cond, parent_nodes, target_node, 0);

    parents_ = parent_nodes;
    std::swap(table_, new_table);
}

// inline Trick
// cpt_tをtemplate化すれば不要
inline void cpt_t::assign_impl(table_type& new_table, condition_t cond, std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node, std::size_t const n) const
{
    if(n >= parent_nodes.size())
    {
        // 挿入してみる(あれば上書きしない)
        new_table.emplace(cond, std::vector<double>(target_node->selectable_num));
        return;
    }
    else
    {
        for(std::size_t i = 0; i < parent_nodes.at(n)->selectable_num; ++i)
        {
            cond[parent_nodes.at(n)] = i;
            assign_impl(new_table, cond, parent_nodes, target_node, n + 1);
        }
    }
}

} // namespace bn

#endif // #ifndef BNI_GRAPH_HPP
