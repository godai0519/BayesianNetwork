#ifndef BNI_LEARNING_STEPWISE_STRUCTURE_HC_HPP
#define BNI_LEARNING_STEPWISE_STRUCTURE_HC_HPP

#include <array>
#include <string>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/hash.hpp>
#include <bayesian/evaluation/transinformation.hpp>
#include <bayesian/learning/utility.hpp>

namespace bn {
namespace learning {

class mutual_information_holder {
public:
    using pair_dictionary_type = std::unordered_map<std::pair<vertex_type, vertex_type>, double, bn::hash<std::pair<vertex_type, vertex_type>>>;

    mutual_information_holder(bn::sampler const& sampling)
        : sampling_(sampling), entropy_machine_(), mutual_machine_()
    {
    }
    
    double calculate_entropy(vertex_type const& node)
    {
        // 存在しているか
        auto it = entropy_dic_.find(node);
        if(it != entropy_dic_.end()) return it->second;

        // 新たに計算して返す
        auto value = entropy_machine_(sampling_, node);
        entropy_dic_.emplace(node, value);
        return value;
    }
    
    double calculate_joint_entropy(vertex_type const& lhs, vertex_type const& rhs)
    {
        auto jointed = make_vertex_pair(lhs, rhs);

        // 存在しているか
        auto it = joint_entropy_dic_.find(jointed);
        if(it != joint_entropy_dic_.end()) return it->second;

        // 新たに計算して返す
        auto value = entropy_machine_(sampling_, {jointed.first, jointed.second});
        joint_entropy_dic_.emplace(std::move(jointed), value);
        return value;
    }

    double calculate_similarity(vertex_type const& lhs, vertex_type const& rhs)
    {
        auto jointed = make_vertex_pair(lhs, rhs);

        // 存在しているか
        auto it = similarity_dic_.find(jointed);
        if(it != similarity_dic_.end()) return it->second;

        // 新たに計算して返す
        auto value = mutual_machine_(
            calculate_entropy(lhs),
            calculate_entropy(rhs),
            calculate_joint_entropy(lhs, rhs)
            );
        similarity_dic_.emplace(std::move(jointed), value);
        return value;
    }

    void delete_entropy(vertex_type const& node)
    {
        entropy_dic_.erase(node);
    }
    
    void delete_joint_entropy(vertex_type const& lhs, vertex_type const& rhs)
    {
        joint_entropy_dic_.erase(make_vertex_pair(lhs, rhs));
    }

    void delete_similarity(vertex_type const& lhs, vertex_type const& rhs)
    {
        similarity_dic_.erase(make_vertex_pair(lhs, rhs));
    }
    
private:
    std::pair<vertex_type, vertex_type> make_vertex_pair(vertex_type const& lhs, vertex_type const& rhs)
    {
        return (lhs < rhs) ? std::make_pair(lhs, rhs)
                           : std::make_pair(rhs, lhs);
    }

    sampler const& sampling_;
    bn::evaluation::entropy const entropy_machine_;
    bn::evaluation::mutual_information const mutual_machine_;
    std::unordered_map<vertex_type, double> entropy_dic_;
    pair_dictionary_type joint_entropy_dic_;
    pair_dictionary_type similarity_dic_;
};

template<class Eval, template<class> class BetweenLearning>
class stepwise_structure_hc {
public:
    using cluster_type = std::vector<vertex_type>;
    using similarity_type = std::tuple<const cluster_type*, const cluster_type*, double>;
    using Similarity = bn::evaluation::mutual_information;

    stepwise_structure_hc(bn::sampler const& sampling)
        : sampling_(sampling), eval_(sampling), learning_machine_(sampling_), engine_(make_engine<std::mt19937>()), info_holder_(sampling)
    {
    }

    double operator()(graph_t& graph)
    {
        // graphの初期化
        graph.erase_all_edge();

        // 初期クラスタと初期類似度を得る
        auto clusters = initial_clustering(graph.vertex_list()); // 初期クラスタ
        auto similarities = initial_similarity(clusters);        // 初期類似度

        // クラスタ間学習(結合)
        auto const score = learning_between_clusters(graph, clusters, similarities);

        return score;
    }

private:
    // 初期クラスタリングを行い，クラスタリング結果を返す
    // 第1引数: クラスタリング対象のノード集合
    std::vector<cluster_type> initial_clustering(std::vector<vertex_type> const& nodes)
    {
        // クラスタ集合
        std::vector<cluster_type> clusters(nodes.size());

        // 1ノード1クラスタ
        for(std::size_t i = 0; i < nodes.size(); ++i)
        {
            clusters[i].push_back(nodes[i]);
        }

        return clusters;
    }

    std::vector<similarity_type> initial_similarity(std::vector<cluster_type> const& clusters)
    {
        std::vector<similarity_type> similarities;
        for(std::size_t i = 0; i < clusters.size(); ++i)
        {
            for(std::size_t j = i + 1; j < clusters.size(); ++j)
            {
                auto const& i_cluster = clusters[i];
                auto const& j_cluster = clusters[j];
                similarities.push_back(make_similarity_tuple(i_cluster, j_cluster));
            }
        }
        return similarities;
    }

    similarity_type make_similarity_tuple(cluster_type const& lhs, cluster_type const& rhs)
    {
        // 2クラスタ間のノードのそれぞれの組み合わせ数
        auto const combination_num = lhs.size() * rhs.size();

        // 類似度計算
        double value = 0;
        for(auto const& lhs_nodes : lhs)
        {
            for(auto const& rhs_nodes : rhs)
            {
                // 数で割って足す(平均)
                value += info_holder_.calculate_similarity(lhs_nodes, rhs_nodes) / combination_num;
            }
        }

        // アドレスが小さいクラスタを先にして返す
        return (&lhs < &rhs) ? std::make_tuple(&lhs, &rhs, value)
                             : std::make_tuple(&rhs, &lhs, value);
    }

    // クラスタ間学習を行う
    double learning_between_clusters(
        graph_t& graph,
        std::vector<cluster_type>& clusters,
        std::vector<similarity_type>& similarities
        )
    {
        double score = std::numeric_limits<double>::max();

        while(clusters.size() != 1)
        {
            //
            // 階層的構造学習
            //

            // 最も似ているクラスタ間
            auto const most_similar = std::max_element(
                similarities.begin(), similarities.end(),
                [](similarity_type const& lhs, similarity_type const& rhs){ return std::get<2>(lhs) < std::get<2>(rhs); }
                );

            // most_similarのどちらが親か
            std::uniform_int_distribution<std::size_t> binary_dist(0, 1);
            const cluster_type *parent, *child;
            if(binary_dist(engine_)) std::tie(parent, child, std::ignore) = *most_similar;
            else                     std::tie(child, parent, std::ignore) = *most_similar;

            // parentとchildのindex
            std::size_t const parent_index =
                std::distance(clusters.begin(), std::find(clusters.begin(), clusters.end(), *parent));
            std::size_t const child_index  =
                std::distance(clusters.begin(), std::find(clusters.begin(), clusters.end(), *child));
            
            // learning
            score = learning_machine_.learn_with_hint(graph, *parent, *child);

            // 合成クラスタ
            cluster_type new_cluster;
            new_cluster.reserve(parent->size() + child->size());
            new_cluster.insert(new_cluster.end(), parent->cbegin(), parent->cend());
            new_cluster.insert(new_cluster.end(), child->cbegin(), child->cend());

            // 前のクラスタに関する類似度を削除
            similarities.erase(
                std::remove_if(
                    similarities.begin(), similarities.end(),
                    [parent, child](similarity_type const& similarity)
                    {
                        // 多分設計ミスってる
                        return
                            std::get<0>(similarity) == parent || std::get<0>(similarity) == child  ||
                            std::get<1>(similarity) == parent || std::get<1>(similarity) == child;
                    }),
                similarities.end()
                );

            // 前のクラスタを消して，挿入
            clusters.erase(clusters.begin() + parent_index);
            clusters.erase(clusters.begin() + child_index + (parent_index < child_index ? -1 : 0));
            clusters.push_back(std::move(new_cluster));

            //
            // 確率的枝刈り (Yet)
            //

            // 合成クラスタと各クラスタ間の類似度を更新
            cluster_type const& inserted_cluster = clusters.back(); // 挿入後のnew_clusterの参照
            std::for_each(
                clusters.begin(), clusters.end() - 1,
                [this, &inserted_cluster, &similarities](cluster_type const& cluster)
                {
                    similarities.push_back(make_similarity_tuple(inserted_cluster, cluster));
                });
        }

        return score;
    }

    double hierarchical_clustering(
        graph_t& graph,
        std::vector<cluster_type>& clusters,
        std::vector<similarity_type>& similarities
        )
    {
        // 最も似ているクラスタ間
        auto const most_similar = std::max_element(
            similarities.begin(), similarities.end(),
            [](similarity_type const& lhs, similarity_type const& rhs){ return std::get<2>(lhs) < std::get<2>(rhs); }
            );

        // most_similarのどちらが親か
        std::uniform_int_distribution<std::size_t> binary_dist(0, 1);
        const cluster_type *parent, *child;
        if(binary_dist(engine_)) std::tie(parent, child, std::ignore) = *most_similar;
        else                     std::tie(child, parent, std::ignore) = *most_similar;

        // parentとchildのindex
        std::size_t const parent_index =
            std::distance(clusters.begin(), std::find(clusters.begin(), clusters.end(), *parent));
        std::size_t const child_index  =
            std::distance(clusters.begin(), std::find(clusters.begin(), clusters.end(), *child));
            
        // learning
        return learning_machine_.learn_with_hint(graph, *parent, *child);
    }

    sampler sampling_;
    Eval eval_;
    BetweenLearning<Eval> learning_machine_;
    std::mt19937 engine_;
    mutual_information_holder info_holder_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_STEPWISE_STRUCTURE_HC_HPP
