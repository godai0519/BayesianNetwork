#ifndef BNI_LEARNING_STEPWISE_STRUCTURE_HPP
#define BNI_LEARNING_STEPWISE_STRUCTURE_HPP

#include <string>
#include <random>
#include <bayesian/graph.hpp>
#include <bayesian/sampler.hpp>
#include <bayesian/learning/utility.hpp>

namespace bn {
namespace learning {

template<class Eval>
class stepwise_structure {
public:
    using cluster_type = std::vector<vertex_type>;

    stepwise_structure(bn::sampler const& sampling)
        : sampling_(sampling), eval_(sampling), engine_(make_engine<std::mt19937>())
    {
    }

    double operator()(graph_t& graph, std::size_t const& initial_cluster_size)
    {
        // graphの初期化
        graph.erase_all_edge();

        // 初期クラスタを得る
        auto clusters = initial_clustering(graph.vertex_list(), initial_cluster_size);

        // 初期クラスタのクラスタ内学習
        learning_intercluster(graph, clusters);

        // クラスタ間学習(結合)
        auto const score = learning_between_clusters(graph, clusters);

        return score;
    }

private:
    // 初期クラスタリングを行い，クラスタリング結果を返す
    // 第1引数: クラスタリング対象のノード集合
    // 第2引数: クラスタリングの最大サイズ
    std::vector<cluster_type> initial_clustering(std::vector<vertex_type> const& nodes, std::size_t const& initial_cluster_size)
    {
        // クラスタ集合のサイズ
        std::size_t const cluster_num =
            (nodes.size() / initial_cluster_size) + (nodes.size() % initial_cluster_size ? 1 : 0);
        std::vector<cluster_type> clusters(cluster_num); // クラスタ集合

        // ノードをシャッフルしておく
        auto random_nodes = nodes;
        std::shuffle(random_nodes.begin(), random_nodes.end(), engine_);

        // 順番に投入
        std::size_t i = 0;
        for(auto it = random_nodes.cbegin(), end = random_nodes.cend(); it != end; ++it)
        {
            clusters[i].push_back(*it);
            i = (i + 1) % cluster_num;
        }

        return clusters;
    }

    // クラスタ内学習を行う
    void learning_intercluster(graph_t& graph, std::vector<cluster_type> const& clusters)
    {
        bn::learning::brute_force<Eval> learning_machine(sampling_);

        for(auto const& cluster : clusters) 
            learning_machine(graph, cluster);
    }

    // クラスタ間学習を行う(とりあえずGreedy)
    double learning_between_clusters(graph_t& graph, std::vector<cluster_type>& clusters)
    {
        bn::learning::greedy<Eval> learning_machine(sampling_);
        double score = std::numeric_limits<double>::max();

        // 2クラスタを無作為に選び，結合
        while(clusters.size() != 1)
        {
            // 親と子を選ぶ
            std::uniform_int_distribution<std::size_t> dist(0, clusters.size() - 1);
            std::size_t const child_index = dist(engine_);
            std::size_t parent_index;
            while((parent_index = dist(engine_)) == child_index);

            auto const& parent = clusters[parent_index];
            auto const& child = clusters[child_index];
            score = learning_machine.learn_with_hint(graph, parent, child);

            // 合成クラスタ
            cluster_type new_cluster;
            new_cluster.reserve(parent.size() + child.size());
            new_cluster.insert(new_cluster.end(), parent.cbegin(), parent.cend());
            new_cluster.insert(new_cluster.end(), child.cbegin(), child.cend());

            // 前のクラスタを消して，挿入
            clusters.erase(clusters.begin() + parent_index);
            clusters.erase(clusters.begin() + child_index + (parent_index < child_index ? -1 : 0));
            clusters.push_back(std::move(new_cluster));
        }

        return score;
    }

    sampler sampling_;
    Eval eval_;
    std::mt19937 engine_;
};

} // namespace learning
} // namespace bn

#endif // BNI_LEARNING_STEPWISE_STRUCTURE_HPP
