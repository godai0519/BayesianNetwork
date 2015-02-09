#include <algorithm>
#include <limits>
#include "bayesian/graph.hpp"

namespace bn {

graph_t& graph_t::operator=(graph_t const& rhs)
{
    graph_t(rhs).swap(*this);
    return *this;
}

#if defined(_MSC_VER) && _MSC_VER < 1900
graph_t::graph_t(graph_t && other)
{
    other.swap(*this);
}

graph_t& graph_t::operator=(graph_t&& rhs)
{
    rhs.swap(*this);
    return *this;
}
#else
graph_t::graph_t(graph_t && other) = default;
graph_t& graph_t::operator=(graph_t&& rhs) = default;
#endif


#if defined(_MSC_VER) && _MSC_VER < 1900
void graph_t::swap(graph_t& other)
#else
void graph_t::swap(graph_t& other) noexcept
#endif
{
    std::swap(vertex_list_, other.vertex_list_);
    std::swap(edge_list_, other.edge_list_);
    std::swap(adjacent_list_, other.adjacent_list_);
}

#if defined(_MSC_VER) && _MSC_VER < 1900
void swap(graph_t& lhs, graph_t& rhs)
#else
void swap(graph_t& lhs, graph_t& rhs) noexcept
#endif
{
    lhs.swap(rhs);
}

std::vector<vertex_type> const& graph_t::vertex_list() const
{
    return vertex_list_;
}

std::vector<edge_type> const& graph_t::edge_list() const
{
    return edge_list_;
}

vertex_type graph_t::add_vertex()
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

edge_type graph_t::add_edge(vertex_type const& from, vertex_type const& to)
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

bool graph_t::erase_vertex(vertex_type const& v)
{
    auto const index = index_search(v);
    if(index >= vertex_list_.size()) return false;

    vertex_list_.erase(vertex_list_.begin() + index);
    adjacent_list_.erase(adjacent_list_.begin() + index);
    for(auto& line : adjacent_list_) line.erase(line.begin() + index);
    return true;
}

bool graph_t::erase_edge(edge_type const& e)
{
    auto index = edge_search(e);
    if(index.first >= vertex_list_.size() || index.second >= vertex_list_.size()) return false;

    edge_list_.erase(std::remove(edge_list_.begin(), edge_list_.end(), e), edge_list_.end());
    adjacent_list_[index.first][index.second] = nullptr;
    return true;
}

bool graph_t::erase_all_vertex()
{
    vertex_list_.clear();
    edge_list_.clear();
    adjacent_list_.clear();
    return true;
}

bool graph_t::erase_all_edge()
{
    edge_list_.clear();
    adjacent_list_.assign(vertex_list_.size(), std::vector<edge_type>(vertex_list_.size(), nullptr));
    return true;
}

bool graph_t::change_edge_direction(edge_type const& e)
{
    auto to   = target(e);
    auto from = source(e);

    if (erase_edge(e) == false)
        return false;

    if (add_edge(to, from) == nullptr)
    {
        //ダメだったら戻す
        add_edge(from, to);
        return false;
    }
    return true;
}

std::vector<edge_type> graph_t::out_edges(vertex_type const& from) const
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

std::vector<vertex_type> graph_t::out_vertexes(vertex_type const& from) const
{
    std::vector<vertex_type> ret;
    for(auto const& edge : out_edges(from))
    {
        ret.push_back(target(edge));
    }
    return ret;
}

std::vector<edge_type> graph_t::in_edges(vertex_type const& to) const
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

std::vector<vertex_type> graph_t::in_vertexes(vertex_type const& to) const
{
    std::vector<vertex_type> ret;
    for(auto const& edge : in_edges(to))
    {
        ret.push_back(source(edge));
    }
    return ret;
}

vertex_type graph_t::source(edge_type const& edge) const
{
    auto const index = edge_search(edge);
    auto const tmp = std::numeric_limits<std::size_t>::max();
    if(index.first == tmp || index.second == tmp) return nullptr;

    return vertex_list_[index.first];
}

vertex_type graph_t::target(edge_type const& edge) const
{
    auto const index = edge_search(edge);
    auto const tmp = std::numeric_limits<std::size_t>::max();
    if(index.first == tmp || index.second == tmp) return nullptr;

    return vertex_list_[index.second];
}

bool graph_t::is_able_trace(vertex_type const& from, vertex_type const& to) const
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

std::size_t graph_t::index_search(vertex_type const& vertex) const
{
    auto const iterator = std::find(vertex_list_.cbegin(), vertex_list_.cend(), vertex);
    return std::distance(vertex_list_.cbegin(), iterator);
}

std::pair<std::size_t, std::size_t> graph_t::edge_search(edge_type const& edge) const
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

} // namespace bn
