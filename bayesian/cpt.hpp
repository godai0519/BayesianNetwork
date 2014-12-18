#ifndef BNI_CPT_HPP
#define BNI_CPT_HPP

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
#else
template<> struct hash<bn::condition_t> : public __hash_base<std::size_t, bn::condition_t> {
#endif
    inline std::size_t operator()(bn::condition_t const& cond) const noexcept
    {
        std::size_t value = 0;
        for(auto const& one : cond)
        {
            hash_combine(value, one.first);
            hash_combine(value, one.second);
        }
        return value;
    }
};

} // namespace std

namespace bn {

class cpt_t {
public:
    typedef std::unordered_map<condition_t, std::vector<double>> table_type;

    explicit cpt_t();
    explicit cpt_t(std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node);

    void assign(std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node);

    // 一部の条件を元にそれに該当するデータを抽出して返す
    table_type filter(condition_t const& cond);

    // 条件として使用できるノードリストを返す
    std::vector<vertex_type> condition_node();

    // 条件が完全一致した確率vectorを返す
    // std::pairのfirstが検索成功したかをのせる
    // firstがtrueのとき，secondには実体への参照が格納される
    // firstがfalseのときのsecondについては未定義
    std::pair<bool, std::vector<double>&> operator[] (condition_t const cond);
    std::pair<bool, std::vector<double> const&> operator[] (condition_t const cond) const;

private:
    void assign_impl(table_type& new_table, condition_t cond, std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node, std::size_t const n) const;

    std::vector<vertex_type> parents_;
    table_type table_;
};

} // namespace bn

#endif // #ifndef BNI_CPT_HPP

