#include "bayesian/graph.hpp"
#include "bayesian/bp.hpp"

namespace bn {

// 辺のlikelihoodを保持し，簡易に前の状態をコピーすることが可能にするため
class likelihood_list {
public:
    class value_type {
    public:
        value_type(vertex_type const& from, vertex_type const& to, edge_type const& edge)
            : from_(from), to_(to), edge_(edge), likelihood_(from->selectable_num, to->selectable_num)
        {
        }
        virtual ~value_type() = default;

        vertex_type const from() const { return from_; }
        vertex_type const to() const { return to_; }
        edge_type const edge() const { return edge_; }

        matrix_type& likelihood()
        {
            return likelihood_;
        }
        matrix_type const& likelihood() const
        {
            return likelihood_;
        }

    private:
        vertex_type from_;
        vertex_type to_;
        edge_type edge_;
        matrix_type likelihood_;
    };

    matrix_type& add_manage(vertex_type const& from, vertex_type const& to, edge_type const& edge)
    {
        data_.emplace_back(from, to, edge);
        return data_.back().likelihood();
    }

    void del_manage(edge_type const& edge)
    {
        auto it = find(edge);
        if(it != data_.end()) data_.erase(it);
    }

    void del_manage(vertex_type const& from, vertex_type const& to)
    {
        auto it = find(from, to);
        if(it != data_.end()) data_.erase(it);
    }

    matrix_type& operator() (edge_type const& edge)
    {
        auto it = find(edge);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (edge_type)");
        return it->likelihood();
    }

    matrix_type const& operator() (edge_type const& edge) const
    {
        auto it = find(edge);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (edge_type) const");
        return it->likelihood();
    }

    matrix_type& operator() (vertex_type const& from, vertex_type const& to)
    {
        auto it = find(from, to);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (vertex_type,vertex_type)");
        return it->likelihood();
    }

    matrix_type const& operator() (vertex_type const& from, vertex_type const& to) const
    {
        auto it = find(from, to);
        if(it == data_.end()) throw std::runtime_error("likelihood_list: operator() (vertex_type,vertex_type) const");
        return it->likelihood();
    }

private:
    std::vector<value_type>::iterator find(edge_type const& edge)
    {
        auto it = data_.begin();
        while(it != data_.end())
        {
            if(it->edge() == edge) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type>::iterator find(vertex_type const& from, vertex_type const& to)
    {
        auto it = data_.begin();
        while(it != data_.end())
        {
            if(std::tie(it->from(), it->to()) == std::tie(from, to)) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type>::const_iterator find(edge_type const& edge) const
    {
        auto it = data_.cbegin();
        while(it != data_.cend())
        {
            if(it->edge() == edge) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type>::const_iterator find(vertex_type const& from, vertex_type const& to) const
    {
        auto it = data_.cbegin();
        while(it != data_.cend())
        {
            if(std::tie(it->from(), it->to()) == std::tie(from, to)) break;
            else ++it;
        }
        return it;
    }

    std::vector<value_type> data_;
};

matrix_type bp::operator()(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    // 前後の要素に伝播させる
    auto const e_minus = propagate_forward(graph, node, condition);
    auto const e_plus = propagate_backward(graph, node, condition);

    // 掛け算
    auto const elem_num = node->selectable_num;
    double sum = 0.0;
    matrix_type mat(elem_num, 1);
    for(std::size_t i = 0; i < e_minus.height(); ++i)
    {
        double const product = e_minus[i][0] * e_plus[0][i];
        sum += product;
        mat[i][0] = product;
    }

    // 正規化
    for(std::size_t i = 0; i < e_minus.height(); ++i)
    {
        mat[i][0] /= sum;
    }

    return mat;
}

std::pair<bool, int> find_condition(
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    for(auto const& c : condition)
    {
        if(c.first == node)
        {
            return std::make_pair(true, c.second);
        }
    }

    return std::make_pair(false, 0);
}

// 下流要素の確率推論
matrix_type bp::propagate_forward(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    auto const elem_num = node->selectable_num;
    matrix_type mat(elem_num, 1, 1);

    // node ∈ condition
    auto is_condition = find_condition(node, condition);
    if(is_condition.first)
    {
        for(int i = 0; i < elem_num; ++i)
        {
            mat[i][0] = (i == is_condition.second) ? 1 : 0;
        }
        return mat;
    }

    // conditionに含まれないから伝播 (e-要素)
    auto const out_edges = graph.out_edges(node);
    if(!out_edges.empty())
    {
        for(auto const& edge : out_edges)
        {
            if(!edge->likelihood.first) throw std::runtime_error("no set edge of likelihood");
            mat = mat % (edge->likelihood.second * propagate_forward(graph, graph.target(edge), condition));
        }

        return mat;
    }

    // 末端は全ての確率が等しいとする
    for(int i = 0; i < elem_num; ++i)
    {
        mat[i][0] = 1.0;
    }
    return mat;
}

// 上流要素の確率推論
matrix_type bp::propagate_backward(
    graph_t const& graph,
    vertex_type const& node,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    auto const elem_num = node->selectable_num;
    matrix_type mat(1, elem_num, 1);

    // node ∈ condition
    auto is_condition = find_condition(node, condition);
    if(is_condition.first)
    {
        for(int i = 0; i < elem_num; ++i)
        {
            mat[0][i] = (i == is_condition.second) ? 1 : 0;
        }
        return mat;
    }

    // conditionに含まれないから伝播 (e+要素)
    auto const in_edges = graph.in_edges(node);
    if(!in_edges.empty())
    {
        for(auto const& edge : in_edges)
        {
            if(!edge->likelihood.first) throw std::runtime_error("no set edge of likelihood");
            mat = mat % (propagate_backward(graph, graph.source(edge), condition) * edge->likelihood.second);
        }

        return mat;
    }

    // 最上位ノードは事前確率を割り当てる
    auto& e = node->evidence;
    if(e.first)
    {
        return e.second;
    }
    else
    {
        throw std::runtime_error("highest node doesn't have prior probability.");
    }
}

} // namespace bn

