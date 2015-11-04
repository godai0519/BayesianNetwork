#ifndef BNI_LEARNING_STEPWISE_STRUCTURE_HC_HPP
#define BNI_LEARNING_STEPWISE_STRUCTURE_HC_HPP

#include <array>
#include <string>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/hash.hpp>
#include <bayesian/evaluation/transinformation.hpp>
#include <bayesian/utility.hpp>

namespace bn {
namespace learning {

template<class Eval, template<class> class BetweenLearning, class PruningProbExpr>
class stepwise_structure_hc {
public:
    using cluster_type = std::shared_ptr<std::vector<vertex_type>>;
    using similarity_type = std::tuple<std::tuple<cluster_type, cluster_type>, double, int>;
    using Similarity = bn::evaluation::mutual_information;

    stepwise_structure_hc(bn::sampler const& sampling)
        : sampling_(sampling), learning_machine_(sampling_), mutual_information_machine_(), engine_(make_engine<std::mt19937>())
    {
    }

    // 階層的クラスタリング及び確率的枝刈りを用いた段階的構造学習法の実行
    // graph: グラフ構造(どんな構造が入っていたとしても，クリアされる)
    // alpha: 枝刈りが実行される確率に関する係数
    double operator()(graph_t& graph, double const alpha)
    {
        // graphの初期化
        graph.erase_all_edge();

        // 初期クラスタと初期類似度を初期化
        initial_clustering(graph.vertex_list()); // 初期クラスタ
        initial_similarities();                  // 初期類似度

        // クラスタ間学習(結合)
        auto const score = learning_between_clusters(graph, alpha);
        return score;
    }

private:
    // 初期クラスタリングを行い，クラスタリング結果をメンバ変数clusters_に格納する
    // 第1引数: クラスタリング対象のノード集合
    void initial_clustering(std::vector<vertex_type> const& nodes)
    {
        // クラスタ集合初期化
        clusters_.clear();
        clusters_.reserve(nodes.size());

        // 1ノード1クラスタ
        for(auto const& node : nodes)
        {
            auto cluster = std::make_shared<cluster_type::element_type>();
            cluster->push_back(node);
            clusters_.push_back(std::move(cluster));
        }
    }

    // 初期類似度計算を行い，類似度をメンバ変数similarities_に格納する
    void initial_similarities()
    {
        // 類似度集合初期化
        similarities_.clear();
        similarities_.reserve(clusters_.size() * clusters_.size());

        // 類似度平均初期化
        average_similar_ = 0.0;
        auto const max_edge_num = clusters_.size() * (clusters_.size() - 1) / 2;

        // 全クラスタペアについて，類似度計算
        for(std::size_t i = 0; i < clusters_.size(); ++i)
        {
            for(std::size_t j = i + 1; j < clusters_.size(); ++j)
            {
                auto const& i_cluster = clusters_[i];
                auto const& j_cluster = clusters_[j];
                assert(i_cluster->size() == 1 && j_cluster->size() == 1); // 初期状態は1ノードであるはずだから

                auto&& similarity = make_similarity_tuple(i_cluster, j_cluster);
                average_similar_ += std::get<1>(similarity) / max_edge_num;
                similarities_.push_back(std::move(similarity));
            }
        }
    }

    // 指定したsimilarityにclusterが関与しているか(clusterに関する類似度か)どうかを返す
    bool is_related(similarity_type const& similarity, cluster_type const& cluster)
    {
        auto const& connection = std::get<0>(similarity);
        return std::get<0>(connection) == cluster || std::get<1>(connection) == cluster;
    }

    // 指定したsimilarityがlhsとrhsに関する類似度かどうかを返す
    bool is_connected(similarity_type const& similarity, cluster_type const& lhs, cluster_type const& rhs)
    {
        auto const& connection = std::get<0>(similarity);
        return std::get<0>(connection) == std::min(lhs, rhs) && std::get<1>(connection) == std::max(lhs, rhs);
    }

    // 引数の2つのクラスタを1つのクラスタに結合する
    // clusters_に追加されていれば，clusters_から削除する(副作用)
    cluster_type combine_clusters(cluster_type const& lhs, cluster_type const& rhs)
    {
        // 合成クラスタ
        auto new_cluster = std::make_shared<cluster_type::element_type>();
        new_cluster->reserve(lhs->size() + rhs->size());
        new_cluster->insert(new_cluster->end(), lhs->cbegin(), lhs->cend());
        new_cluster->insert(new_cluster->end(), rhs->cbegin(), rhs->cend());

        // 前のクラスタを消す
        clusters_.erase(std::find(clusters_.begin(), clusters_.end(), lhs));
        clusters_.erase(std::find(clusters_.begin(), clusters_.end(), rhs));

        return new_cluster;
    }

    // メンバ変数similarities_から最も順位の高いものを取り出し，無作為に親子を決定する
    similarity_type most_similarity()
    {
        // 最も似ているクラスタ間
        auto const most_similar = std::max_element(
            similarities_.begin(), similarities_.end(),
            [](similarity_type const& lhs, similarity_type const& rhs)
            {
                return std::get<2>(lhs) < std::get<2>(rhs)
                    || !(std::get<2>(rhs) < std::get<2>(lhs)) && std::get<1>(lhs) < std::get<1>(rhs);
            });

        // コピーして元のクラスタ間を消す
        auto result = *most_similar;
        similarities_.erase(most_similar);

        // most_similarのどちらが親か(一定条件で入れ替え)
        std::uniform_int_distribution<std::size_t> binary_dist(0, 1);
        auto pair = std::get<0>(result);
        if(binary_dist(engine_)) std::swap(std::get<0>(pair), std::get<1>(pair));

        return result;
    }
    
    std::tuple<cluster_type, cluster_type> make_cluster_tuple(cluster_type const& lhs, cluster_type const& rhs)
    {
        // アドレスが小さいクラスタを先にして返す
        return (lhs < rhs) ? std::make_tuple(lhs, rhs)
                           : std::make_tuple(rhs, lhs);
    }

    // 与えられた2つのクラスタ間の類似度を求め，similarity_typeとして返す
    similarity_type make_similarity_tuple(cluster_type const& lhs, cluster_type const& rhs)
    {
        // 2クラスタ間のノードのそれぞれの組み合わせ数
        auto const combination_num = lhs->size() * rhs->size();

        // 類似度計算
        double value = 0;
        for(auto const& lhs_nodes : *lhs)
        {
            for(auto const& rhs_nodes : *rhs)
            {
                // 数で割って足す(平均)
                value += mutual_information_machine_(sampling_, lhs_nodes, rhs_nodes) / combination_num;
            }
        }

        return std::make_tuple(make_cluster_tuple(lhs, rhs), value, 1);
    }

    // クラスタ間学習を行う
    // graph, alpha: operator()参照
    double learning_between_clusters(graph_t& graph, double const alpha)
    {
        double score = std::numeric_limits<double>::max();

        while(clusters_.size() != 1 && !similarities_.empty())
        {
            //
            // 階層的構造学習 部分
            //

            // 結合対象を得る
            auto const similarity_target = most_similarity();
            if(std::get<2>(similarity_target) != 1) break;
            
            // learning
            auto const combine_target = std::get<0>(similarity_target);
            auto const parent = std::get<0>(combine_target);
            auto const child  = std::get<1>(combine_target);
            score = learning_machine_.learn_with_hint(graph, *parent, *child);

            // クラスタ合成
            auto combined_cluster = combine_clusters(parent, child);

            //
            // 確率的枝刈り 部分
            //
            stochastic_pruning(alpha, combined_cluster, similarity_target);
            clusters_.push_back(combined_cluster);
        }

        return score;
    }

    // 確率的枝刈りを行う
    // alpha: operator()に準ずる
    // new_cluster: 結合後のクラスタを示す
    // old_connection: 結合前の2クラスタ間のsimilarity_typeを示す
    void stochastic_pruning(
        double const alpha,
        cluster_type const& new_cluster, similarity_type const& old_similarity
        )
    {
        // Initialize
        std::unordered_map<cluster_type, std::tuple<double, std::vector<similarity_type>>> next_similarity;
        for(auto const& cluster : clusters_)
            next_similarity[cluster] = std::make_tuple(0.0, std::vector<similarity_type>());

        // 旧クラスタ類似度より新クラスタ類似度を算出し，消去する
        for(auto it = similarities_.begin(); it != similarities_.end();)
        {
            auto const& connection = std::get<0>(*it);
            auto const& old_connection = std::get<0>(old_similarity);
            if(std::get<0>(connection) == std::get<0>(old_connection) || std::get<0>(connection) == std::get<1>(old_connection))
            {
                auto const& old_cluster = std::get<0>(connection);
                auto const& target_cluster = std::get<1>(connection);
                auto const similarity = std::get<1>(*it);
                
                // 格納と削除
                auto& container = next_similarity[target_cluster];
                std::get<0>(container) += similarity * old_cluster->size() / new_cluster->size();
                std::get<1>(container).push_back(*it);
                it = similarities_.erase(it);
            }
            else if(std::get<1>(connection) == std::get<0>(old_connection) || std::get<1>(connection) == std::get<1>(old_connection))
            {
                auto const& old_cluster = std::get<1>(connection);
                auto const& target_cluster = std::get<0>(connection);
                auto const similarity = std::get<1>(*it);
                
                // 格納と削除
                auto& container = next_similarity[target_cluster];
                std::get<0>(container) += similarity * old_cluster->size() / new_cluster->size();
                std::get<1>(container).push_back(*it);
                it = similarities_.erase(it);
            }
            else ++it;
        }

        // 新クラスタ類似度より枝刈りを実行する
        for(std::size_t i = 0, end = clusters_.size(); i < end; ++i)
        {
            auto const& prunning_target = next_similarity[clusters_[i]];

            unsigned char flag = 0;
            flag |= std::get<2>(std::get<1>(prunning_target)[0]);
            flag |= std::get<2>(std::get<1>(prunning_target)[1]) * 2;

            double const prunning_probability =
                (flag == 3) ? pruning_probability_.p1(alpha, average_similar_, std::get<1>(old_similarity), std::get<0>(prunning_target), std::get<1>(std::get<1>(prunning_target)[0]), std::get<1>(std::get<1>(prunning_target)[1])) :
                (flag == 2) ? pruning_probability_.p2(alpha, average_similar_, std::get<1>(old_similarity), std::get<0>(prunning_target), std::get<1>(std::get<1>(prunning_target)[1]), std::get<1>(std::get<1>(prunning_target)[0])) :
                (flag == 1) ? pruning_probability_.p2(alpha, average_similar_, std::get<1>(old_similarity), std::get<0>(prunning_target), std::get<1>(std::get<1>(prunning_target)[0]), std::get<1>(std::get<1>(prunning_target)[1])) :
                              pruning_probability_.p3(alpha, average_similar_, std::get<1>(old_similarity), std::get<0>(prunning_target), std::get<1>(std::get<1>(prunning_target)[0]), std::get<1>(std::get<1>(prunning_target)[1]));

            // 確率により，probabilityの確率で枝刈り
            std::uniform_real_distribution<double> probability_dist(0.0, 1.0);
            int const is_connection = (probability_dist(engine_) < prunning_probability) ? 0 : 1; // 枝刈りed = 0
            
            similarities_.emplace_back(make_cluster_tuple(clusters_[i], new_cluster), std::get<0>(prunning_target), is_connection);
        }
    }

    sampler const& sampling_;
    BetweenLearning<Eval> learning_machine_;
    bn::evaluation::mutual_information mutual_information_machine_;
    std::mt19937 engine_;

    PruningProbExpr pruning_probability_;

    std::vector<cluster_type> clusters_;
    std::vector<similarity_type> similarities_;
    double average_similar_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_STEPWISE_STRUCTURE_HC_HPP
