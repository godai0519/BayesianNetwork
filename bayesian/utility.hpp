#ifndef BNI_UTILITY_HPP
#define BNI_UTILITY_HPP

#include <string>
#include <vector>
#include <bayesian/graph.hpp>

namespace bn {
namespace io {

template<class InputStream>
std::vector<std::string> stream_to_lines(InputStream& is)
{
    std::vector<std::string> result;
    std::string line;
    while(std::getline(is, line))
    {
        if(!line.empty() && line.back() == '¥r') line.pop_back();
        result.push_back(line);
    }

    return result;
}

// template使って書きたさしかない
// つーかゴリ押し
class dsc {
public:
    graph_t parse(std::vector<std::string> data);
    graph_t from_file(std::string const& filename);
    graph_t from_data(std::string const& data);

private:
    typedef std::vector<std::string>::iterator LineIterator;

    void parse_header(LineIterator& it, LineIterator const& end);
    void parse_node(LineIterator& it, LineIterator const& end);
    void parse_cpt(LineIterator& it, LineIterator const& end);

    graph_t graph_;
    std::unordered_map<std::string, vertex_type> dictionary_;
};


} // namespace io
} // namespace bn

#endif // #ifndef BNI_UTILITY_HPP
