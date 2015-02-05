#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <bayesian/utility.hpp>

namespace bn {
namespace io {

graph_t dsc::parse(std::vector<std::string> data)
{
    // 空行の削除及び行頭の空白を削除
    data.erase(std::remove_if(data.begin(), data.end(), [](std::string const& str){ return str.empty(); }), data.end());
    for(auto& str : data) // boost::algorithm::trim
    {
        auto it = str.begin();
        while(it != str.end() && (!std::isprint(*it) || std::isspace(*it))) ++it;
        str.erase(str.begin(), it);
    }

    // Parse Start
    graph_.erase_all_vertex();
    dictionary_.clear();

    auto it = data.begin();
    auto const end = data.end();
    while(it != end)
    {
        // まだswitchパターン使ったほうがまともだわ
        if(it->substr(0, 14) == "belief network")
        {
            parse_header(it, end);
        }
        else if(it->substr(0, 4) == "node")
        {
            parse_node(it, end);
        }
        else if(it->substr(0, 11) == "probability")
        {
            parse_cpt(it, end);
        }
        else ++it;
    }

    return graph_;
}

graph_t dsc::from_file(std::string const& filename)
{
    std::ifstream ifs(filename);
    auto raw_lines = stream_to_lines(ifs);
    return parse(raw_lines);
}

graph_t dsc::from_data(std::string const& data)
{
    std::istringstream iss(data);
    auto raw_lines = stream_to_lines(iss);
    return parse(raw_lines);
}

void dsc::parse_header(LineIterator& it, LineIterator const& end)
{
    assert(it->substr(0, 14) == "belief network");

    std::string key = it->substr(16); // key"
    key.pop_back();  // key

    // TODO: Graph name(key)の使い方

    ++it;
}

void dsc::parse_node(LineIterator& it, LineIterator const& end)
{
    // 1行目
    assert(it->substr(0, 4) == "node");
    std::string const node_name = it->substr(5);
    dictionary_[node_name] = graph_.add_vertex();
    ++it;

    // 2行目({)
    ++it;

    // 3行目以降
    while(*it != "}")
    {
        if(it->substr(0, 2) == "//")
        {
            ++it;
            continue;
        }

        std::string const key(it->cbegin(), std::find(it->cbegin(), it->cend(), ':'));
        if(key == "type")
        {
            std::string num = it->substr(15);
            num.erase(num.find(']'));
            dictionary_[node_name]->selectable_num = std::stoi(num);
        }
        else
        {
            // Not Implement
        }
        ++it;
    }

    ++it;
}

std::vector<std::string> spirit(std::string const& str, std::string const& delim)
{
    std::vector<std::string> res;
    for(std::string::size_type pos = 0; true;)
    {
        auto const next = str.find(delim, pos);
        res.push_back(str.substr(pos, next - pos));

        if(next == std::string::npos) break;
        pos = next + delim.size();
    }
    return res;
}

void dsc::parse_cpt(LineIterator& it, LineIterator const& end)
{
    // 1行目
    assert(it->substr(0, 11) == "probability");

    std::string condition_str = it->substr(12);
    condition_str.pop_back();

    std::string const node_str = condition_str.substr(0, condition_str.find(" | "));
    std::string const cond_str = (condition_str.size() == node_str.size()) ? "" : condition_str.substr(node_str.size() + 3);
    
    vertex_type const target_node = dictionary_[node_str];
    std::vector<vertex_type> cond_nodes;
    for(auto const& str : spirit(cond_str, ", "))
    {
        auto parent = dictionary_[str];
        if(parent == nullptr) continue;

        cond_nodes.push_back(parent);
        graph_.add_edge(parent, target_node);
    }

    target_node->cpt.assign(cond_nodes, target_node);
    ++it;

    // 2行目
    ++it;
    
    // 3行目以降
    condition_t cond;
    if(cond_nodes.empty())
    {
        while(it->substr(0, 2) == "//") ++it;

        std::vector<double> data;
        for(auto const& str : spirit(*it, ", "))
        {
            data.push_back(std::stod(str));
        }
        
        target_node->cpt[cond].second = data;
        ++it;
    }
    else
    {
        while(*it != "}")
        {
            if(it->substr(0, 2) == "//")
            {
                ++it;
                continue;
            }

            auto line = spirit(*it, ": ");
            line[0].erase(0, 1);
            line[0].pop_back();
            line[1].pop_back();

            std::vector<double> data;
            for(auto const& str : spirit(line[1], ", "))
            {
                data.push_back(std::stod(str));
            }

            auto const cond_parsed = spirit(line[0], ", ");
            for(std::size_t i = 0; i < cond_nodes.size(); ++i)
            {
                cond[cond_nodes[i]] = std::stoi(cond_parsed[i]);
            }

            target_node->cpt[cond].second = data;
            ++it;
        }
    }

    ++it;
}

} // namespace io
} // namespace bn
