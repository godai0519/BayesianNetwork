#ifndef BNI_SERIALIZER_DOT_HPP
#define BNI_SERIALIZER_DOT_HPP

#include <string>
#include <vector>
#include <bayesian/graph.hpp>

namespace bn {
namespace serializer {

class dot {
public:
    template<class OutputStream>
    OutputStream& write(OutputStream& ost, bn::graph_t const& graph, bn::database_t const& data)
    {
        // 冒頭
        ost << "digraph " << data.graph_name << "{\n";

        // ノードを書き出す
        for(auto const& node : graph.vertex_list())
            ost << "    " << get_node_identify(node) << " [label=\"" << data.node_name.at(node->id) << "\"];\n";

        // エッジを書き出す
        for(auto const& edge : graph.edge_list())
        {
            auto const source = graph.source(edge);
            auto const target = graph.target(edge);
            ost << "    " << get_node_identify(source) << " -> " << get_node_identify(target) << ";\n";
        }

        ost << "}";
        return ost;
    }

private:
    inline std::string get_node_identify(bn::vertex_type const& node)
    {
        return "Node" + std::to_string(node->id);
    }
};

} // namespace serializer
} // namespace bn

#endif // #ifndef BNI_SERIALIZER_DOT_HPP

