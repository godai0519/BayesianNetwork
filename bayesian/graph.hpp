#ifndef BNI_GRAPH_HPP
#define BNI_GRAPH_HPP

#include <vector>
#include <cstdint>
#include <memory>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include "matrix.hpp"
#include "cpt.hpp"

namespace bn {

struct vertex_t;
struct edge_t;
typedef std::shared_ptr<vertex_t> vertex_type;
typedef std::shared_ptr<edge_t>   edge_type;

struct vertex_t {
    int id;
    std::size_t selectable_num = 0; // 取りうる値の数
    std::pair<bool, matrix_type> evidence;

    cpt_t cpt;
};

struct edge_t {
};

class graph_t {
public:
    graph_t() = default;
    virtual ~graph_t() = default;

    // copy/move ctor
    graph_t(graph_t const& other) = default;
    graph_t(graph_t && other) = default;

    graph_t& operator=(graph_t const& rhs)
    {
        graph_t(rhs).swap(*this);
        return *this;
    }

    graph_t& operator=(graph_t&& rhs) = default;

    void swap(graph_t& other) noexcept
    {
        std::swap(vertex_list_, other.vertex_list_);
        std::swap(edge_list_, other.edge_list_);
        std::swap(adjacent_list_, other.adjacent_list_);
    }

    friend void swap(graph_t& lhs, graph_t& rhs) noexcept
    {
        lhs.swap(rhs);
    }

    std::vector<vertex_type> const& vertex_list() const
    {
        return vertex_list_;
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
    // 既に辺が存在していた場合，nullptrなshared_ptrを即座に返却する
    // 頂点が存在していなかった場合も同様
    edge_type add_edge(vertex_type const& from, vertex_type const& to)
    {
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

        adjacent_list_[index.first][index.second] = nullptr;
        return true;
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

private:
    std::size_t index_search(vertex_type const& vertex) const
    {
        auto const iterator = std::find(vertex_list_.cbegin(), vertex_list_.cend(), vertex);
        return std::distance(vertex_list_.cbegin(), iterator);
    }

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

} // namespace bn

#endif // #ifndef BNI_MATRIX_HPP

