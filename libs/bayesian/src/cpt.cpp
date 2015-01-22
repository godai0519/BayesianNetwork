#include <algorithm>
#include <iterator>
#include "bayesian/graph.hpp"
#include "bayesian/cpt.hpp"

namespace bn {

cpt_t::cpt_t()
{
}

cpt_t::cpt_t(std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node)
{
    assign(parent_nodes, target_node);
}

void cpt_t::assign(std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node)
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

void cpt_t::assign_impl(table_type& new_table, condition_t cond, std::vector<vertex_type> const& parent_nodes, vertex_type const& target_node, std::size_t const n) const
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

auto cpt_t::filter(condition_t const& cond) -> table_type
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

std::vector<vertex_type> cpt_t::condition_node()
{
    return parents_;
}

std::pair<bool, std::vector<double>&> cpt_t::operator[] (condition_t const cond)
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

std::pair<bool, std::vector<double> const&> cpt_t::operator[] (condition_t const cond) const
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



} // namespace bn

