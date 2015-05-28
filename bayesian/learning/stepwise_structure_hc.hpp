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

// (同時)エントロピーを複数回計算する無駄を省くための辞書ホルダ
// 相互情報量規準
class mutual_information_holder {
public:
    using pair_dictionary_type = std::unordered_map<std::pair<vertex_type, vertex_type>, double, bn::hash<std::pair<vertex_type, vertex_type>>>;

    mutual_information_holder(bn::sampler const& sampling)
        : sampling_(sampling), entropy_machine_(), mutual_machine_()
    {
    }

    // エントロピーを計算して返す
    // node: 調べるノード
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

    // 同時エントロピーを計算して返す
    // lhs, rhs: 調べるノード(順不同)
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

    // 類似度を計算して返す
    // 必要ならエントロピーを内部で計算する
    // lhs, rhs: 調べるノード(順不同)
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

    // 登録済みのエントロピーを消す
    void delete_entropy(vertex_type const& node)
    {
        entropy_dic_.erase(node);
    }

    // 登録済みの同時エントロピーを消す
    void delete_joint_entropy(vertex_type const& lhs, vertex_type const& rhs)
    {
        joint_entropy_dic_.erase(make_vertex_pair(lhs, rhs));
    }

    // 登録済みの類似度を消す
    void delete_similarity(vertex_type const& lhs, vertex_type const& rhs)
    {
        similarity_dic_.erase(make_vertex_pair(lhs, rhs));
    }

private:
    // pairの前半がアドレスの小さいノードになるように正規化する
    // (順不同性を保証する)
    std::pair<vertex_type, vertex_type> make_vertex_pair(vertex_type const& lhs, vertex_type const& rhs)
    {
        return (lhs < rhs) ? std::make_pair(lhs, rhs)
                           : std::make_pair(rhs, lhs);
    }

    // 学習器
    sampler const& sampling_;
    bn::evaluation::entropy const entropy_machine_;
    bn::evaluation::mutual_information const mutual_machine_;

    // 辞書
    std::unordered_map<vertex_type, double> entropy_dic_;
    pair_dictionary_type joint_entropy_dic_;
    pair_dictionary_type similarity_dic_;
};

template<class Eval, template<class> class BetweenLearning>
class stepwise_structure_hc {
public:
    using cluster_type = std::shared_ptr<std::vector<vertex_type>>;
    using similarity_type = std::tuple<cluster_type, cluster_type, double>;
    using Similarity = bn::evaluation::mutual_information;

    stepwise_structure_hc(bn::sampler const& sampling)
        : sampling_(sampling), eval_(sampling), learning_machine_(sampling_), engine_(make_engine<std::mt19937>()), info_holder_(sampling)
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
                average_similar_ += std::get<2>(similarity) / max_edge_num;
                similarities_.push_back(std::move(similarity));
            }
        }
    }

    // 指定したsimilarityにclusterが関与しているか(clusterに関する類似度か)どうかを返す
    bool is_related(similarity_type const& similarity, cluster_type const& cluster)
    {
        return std::get<0>(similarity) == cluster || std::get<1>(similarity) == cluster;
    }

    // 指定したsimilarityがlhsとrhsに関する類似度かどうかを返す
    bool is_connected(similarity_type const& similarity, cluster_type const& lhs, cluster_type const& rhs)
    {
        return std::get<0>(similarity) == std::min(lhs, rhs) && std::get<1>(similarity) == std::max(lhs, rhs);
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
    std::tuple<cluster_type, cluster_type, double> most_similarity()
    {
        // 最も似ているクラスタ間
        auto const most_similar = std::max_element(
            similarities_.begin(), similarities_.end(),
            [](similarity_type const& lhs, similarity_type const& rhs){ return std::get<2>(lhs) < std::get<2>(rhs); }
            );

        // コピーして元のクラスタ間を消す
        auto result = *most_similar;
        similarities_.erase(most_similar);

        // most_similarのどちらが親か(一定条件で入れ替え)
        std::uniform_int_distribution<std::size_t> binary_dist(0, 1);
        if(binary_dist(engine_)) std::swap(std::get<0>(result), std::get<1>(result));

        return result;
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
                value += info_holder_.calculate_similarity(lhs_nodes, rhs_nodes) / combination_num;
            }
        }

        // アドレスが小さいクラスタを先にして返す
        return (lhs < rhs) ? std::make_tuple(lhs, rhs, value)
                           : std::make_tuple(rhs, lhs, value);
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
            similarity_type const combine_target = most_similarity();
            auto const parent = std::get<0>(combine_target);
            auto const child  = std::get<1>(combine_target);

            // learning
            score = learning_machine_.learn_with_hint(graph, *parent, *child);

            // クラスタ合成
            clusters_.push_back(combine_clusters(parent, child));
            auto const& inserted_cluster = clusters_.back();      // 挿入後のnew_clusterの参照


            //
            // 確率的枝刈り 部分
            //
            stochastic_pruning(alpha, inserted_cluster, combine_target);
        }

        return score;
    }

    // 確率的枝刈りを行う
    // alpha: operator()に準ずる
    // new_cluster: 結合後のクラスタを示す
    // old_connection: 結合前の2クラスタ間のsimilarity_typeを示す
    void stochastic_pruning(
        double const alpha,
        cluster_type const& new_cluster, similarity_type const& old_connection
        )
    {
        // 全クラスタと検索
        for(auto const cluster : clusters_)
        {
            // 自分自身及び結合後のクラスタならばパスする
            if(cluster == new_cluster) continue;

            // 探索中クラスタが元のクラスタと幾つ接続されていたか調べる
            std::vector<similarity_type> connection;
            for(auto it = similarities_.begin(); it != similarities_.end(); )
            {
                if(is_connected(*it, cluster, std::get<0>(old_connection)) || is_connected(*it, cluster, std::get<1>(old_connection)))
                {
                    connection.push_back(*it);
                    it = similarities_.erase(it);
                }
                else ++it;
            }

            // 類似度を計算しておく
            auto new_similarity = make_similarity_tuple(new_cluster, cluster);

            // 枝刈り条件
            double probability;
            if(connection.size() == 2)
            {
                probability = std::pow(alpha, std::get<2>(new_similarity) / average_similar_);
            }
            else if(connection.size() == 1)
            {
                probability = std::pow(alpha, std::get<2>(old_connection) / std::get<2>(connection[0]));
            }
            else if(connection.size() == 0)
            {
                // 張らない
                continue;
            }
            else throw std::runtime_error("too connection");

            // 確率により，probabilityの確率で枝刈り
            std::uniform_real_distribution<double> probability_dist(0.0, 1.0);
            if(probability_dist(engine_) < probability) continue; // 枝刈り

            // 刈らない
            similarities_.push_back(std::move(new_similarity));
        }

        // 消えたノードに関する類似度をすべて消す(ゴリ押し)
        for(auto it = similarities_.begin(); it != similarities_.end(); )
        {
            // 消えたノードに関係しているかどうか
            if(is_related(*it, std::get<0>(old_connection)) || is_related(*it, std::get<1>(old_connection)))
            {
                it = similarities_.erase(it);
            }
            else ++it;
        }
    }

    sampler sampling_;
    Eval eval_;
    BetweenLearning<Eval> learning_machine_;
    std::mt19937 engine_;
    mutual_information_holder info_holder_;

    std::vector<cluster_type> clusters_;
    std::vector<similarity_type> similarities_;
    double average_similar_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_STEPWISE_STRUCTURE_HC_HPP
