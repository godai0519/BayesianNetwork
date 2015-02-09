#ifndef BNI_GRAPH_HPP
#define BNI_GRAPH_HPP

#include <vector>
#include <cstdint>
#include <memory>
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
    cpt_t cpt;
};

struct edge_t {}; // tag

class graph_t {
public:
    graph_t() = default;
    virtual ~graph_t() = default;

    // copy/move ctor
    graph_t(graph_t const& other) = default;
    graph_t(graph_t && other);

    graph_t& operator=(graph_t const& rhs);
    graph_t& operator=(graph_t&& rhs);

#if defined(_MSC_VER) && _MSC_VER < 1900
    void swap(graph_t& other);
    friend void swap(graph_t& lhs, graph_t& rhs);
#else
    void swap(graph_t& other) noexcept;
    friend void swap(graph_t& lhs, graph_t& rhs) noexcept;
#endif

    std::vector<vertex_type> const& vertex_list() const;
    std::vector<edge_type> const& edge_list() const;

    // 頂点を生成し，そのshared_ptrを返す
    // 必ず成功する
    vertex_type add_vertex();

    // 引数を元に辺を生成し，そのshared_ptrを返す
    // 既に辺が存在していた場合やDAGが成立しなくなる場合は，nullptrなshared_ptrを即座に返却する
    // 頂点が存在していなかった場合も同様
    edge_type add_edge(vertex_type const& from, vertex_type const& to);

    // 引数を元に一致する頂点を削除する
    // その際，その頂点に入出辺も同時に削除する
    // 削除に成功した場合はtrueを返却する
    bool erase_vertex(vertex_type const& v);

    // 引数を元に一致する辺を削除する
    // 削除に成功した場合はtrueを返却する
    bool erase_edge(edge_type const& e);

    // 全ての頂点を削除する(付随的にerase_all_edgeも起こる)
    bool erase_all_vertex();

    // 全ての辺を削除する
    bool erase_all_edge();

    // 引数で指定された辺を逆向きにする
    bool change_edge_direction(edge_type const& e);

    // 引数の頂点から出て行く辺を列挙する
    std::vector<edge_type> out_edges(vertex_type const& from) const;

    // 引数の頂点から出て行く隣接頂点を列挙する
    std::vector<vertex_type> out_vertexs(vertex_type const& from) const;

    // 引数の頂点へ入っていく辺を列挙する
    std::vector<edge_type> in_edges(vertex_type const& to) const;

    // 引数の頂点へ入っていく隣接頂点を列挙する
    std::vector<vertex_type> in_vertexs(vertex_type const& to) const;

    // 辺の先を探す
    vertex_type source(edge_type const& edge) const;

    // 辺の元を探す
    vertex_type target(edge_type const& edge) const;

    // fromの頂点からtoの頂点が辿れるならばtrue，それ以外ならfalseを返す
    // graphはDAGであることを前提とする
    bool is_able_trace(vertex_type const& from, vertex_type const& to) const;

private:
    // vertex_typeから，vertex_list_内のindexを見つける
    std::size_t index_search(vertex_type const& vertex) const;

    // edge_typeからadjacent_list_内の座標を見つける
    // 該当しないときは
    // {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()}
    std::pair<std::size_t, std::size_t> edge_search(edge_type const& edge) const;

    std::vector<vertex_type> vertex_list_;
    std::vector<edge_type>   edge_list_;
    std::vector<std::vector<edge_type>> adjacent_list_;
};

} // namespace bn

#endif // #ifndef BNI_MATRIX_HPP

