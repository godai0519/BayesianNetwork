#ifndef BNI_SERIALIZER_CSV_HPP
#define BNI_SERIALIZER_CSV_HPP

#include <string>
#include <vector>
#include <bayesian/graph.hpp>
#include <boost/algorithm/string.hpp>

namespace bn {
namespace serializer {

class csv {
public:
    template<class InputStream>
    InputStream& load(InputStream& ist, bn::graph_t& graph)
    {
	    auto const& node_list = graph.vertex_list();

	    std::string line;
	    std::vector<std::string> splited_line;

	    // 1行目は読み飛ばし
	    std::getline(ist, line);

	    // 2行目以降
	    for(int i = 0; i < node_list.size(); ++i)
	    {
		    std::getline(ist, line);
            boost::algorithm::split(splited_line, line, boost::is_any_of(","));
		    for(int j = 0; j < node_list.size(); ++j)
		    {
			    // i と jをつなぐ．見るところはiとj+1
                if(splited_line[j + 1] == "*")
                    graph.add_edge(node_list[i], node_list[j]);
		    }
	    }

        return ist;
    }

    template<class OutputStream>
    OutputStream& write(OutputStream& ost, bn::graph_t const& graph)
    {
        auto const& vertex_list = graph.vertex_list();
        auto const vertex_num = vertex_list.size();

        for(auto const& vertex : vertex_list) ost << "," << vertex->id;
        ost << "\n";

        for(auto const& looking_node : vertex_list)
        {
            ost << looking_node->id;

            auto const& child_nodes = graph.out_vertexes(looking_node);
            for(auto const& candidate_node : vertex_list)
            {
                auto const it = std::find(child_nodes.cbegin(), child_nodes.cend(), candidate_node);
                if(it == child_nodes.cend()) ost << ", ";
                else                         ost << ",*";
            }
            ost << "\n";
        }

        return ost;
    }
};

} // namespace serializer
} // namespace bn

#endif // #ifndef BNI_SERIALIZER_CSV_HPP

