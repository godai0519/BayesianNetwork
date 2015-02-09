#include <algorithm>
#include "bayesian/graph.hpp"
#include "bayesian/rejection_sampling.hpp"

namespace bn {

rejection_sampling::rejection_sampling(graph_t const& graph)
    : graph_(graph)
{
}

rejection_sampling::return_type rejection_sampling::operator()(
    std::vector<std::pair<vertex_type, int>> const& condition,
    int const generate_sample_num
    )
{
    // パターン作る(generate_pattern)
    auto const generated_patterns = generate_pattern(generate_sample_num, condition);

    // 数え上げを行う
    rejection_sampling::return_type result;
    for(auto const node : graph_.vertex_list())
    {
        bn::matrix_type mat(1, node->selectable_num, 0.0);
	    for(auto const& pattern : generated_patterns)
	    {
		    mat[0][pattern.at(node)] += 1.0;
	    }

	    // 全要素をパターン数で割る
	    for(std::size_t i = 0; i < node->selectable_num; ++i)
	    {
		    mat[0][i] /= generated_patterns.size();
	    }

        result[node] = std::move(mat);
    }

    return result;
}

rejection_sampling::pattern_list rejection_sampling::generate_pattern(
    int const num,
    std::vector<std::pair<vertex_type, int>> const& condition
    )
{
    // 条件探査用関数
    auto const is_condition = [&condition](std::unordered_map<vertex_type, int> const& pattern) -> bool
        {
            for(auto const& forcus_cond : condition)
            {
                // ノードを探す．見つからなければfalse
                auto position = pattern.find(forcus_cond.first);
                if(position == pattern.cend()) return false;

                // 条件と一致しなければfalse
                if(position->second != forcus_cond.second) return false;
            }
            return true;
        };

    // 最初の要素をseed_node_listから削除したうえで再帰に移譲(当然ノードがなければ失敗)
    auto const seed_node_list = graph_.vertex_list();
    if(seed_node_list.empty()) return {};


    // num回ループで生成する
    pattern_list result;
    while(result.size() < static_cast<std::size_t>(num))
    {
        auto node_list = seed_node_list;
		std::unordered_map<vertex_type, int> pattern;

        while(!node_list.empty())
        {
            auto const first_node = node_list.back();
            node_list.pop_back();

            choice_pattern(first_node, node_list, pattern);
        }

		if(is_condition(pattern))
		{
			result.push_back(std::move(pattern));
		}
    }

    return result;
}

void rejection_sampling::choice_pattern(
    vertex_type const& current,
    std::vector<vertex_type>& remain_node,
    std::unordered_map<vertex_type, int>& pattern
    )
{
    // 上位ノードが決定していることを確認する．
    // 決定していなかった場合は再帰的処理
    std::unordered_map<vertex_type, int> parent_condition;
    for(auto const& parent : graph_.in_vertexs(current))
    {
        auto const position = std::find(remain_node.begin(), remain_node.end(), parent);
        if(position != remain_node.end())
        {
            // 抜き出して消す
            auto const next_current = *position;
            remain_node.erase(position);

            // 再帰
            choice_pattern(next_current, remain_node, pattern);
        }

        parent_condition[parent] = pattern[parent];
    }

    // CPTから必要データを取得
    auto const optional_cpt_data = current->cpt[parent_condition];
    std::vector<double> condition_probability;
    if(optional_cpt_data.first)
    {
        condition_probability = optional_cpt_data.second;
    }
    else
    {
        throw std::runtime_error("Invalid data from CPT(cannot search data)");
    }

    // メルセンヌ・ツイスタで自身のノード値を生成
    auto const generated_value = probability_generator_();
    double total = 0.0;
    for(int i = 0; current->selectable_num; ++i)
    {
        auto const old_total = total;
        total += condition_probability.at(i);
        if(old_total <= generated_value && generated_value < total)
        {
            pattern[current] = i;
            break;
        }
    }

    return;
}

rejection_sampling::probability_generator::probability_generator()
    : distribution_(0.0, 1.0)
{
    std::random_device rand_dev;
    std::vector<std::uint_least32_t> vec(10);
    std::generate(vec.begin(), vec.end(), std::ref(rand_dev));
    std::seed_seq seed(vec.begin(), vec.end());
    engine_.reset(new std::mt19937(seed));
}

double rejection_sampling::probability_generator::operator() ()
{
    return distribution_(*engine_);
}

} // namespace bn
