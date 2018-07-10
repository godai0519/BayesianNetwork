/**
* @file rejection_sampling.hpp
* @brief Infering (Reasoning) any probability $P(Q|E)$ by using Rejection Sampling.
* @author godai_0519
* @date 02/21/2018
*/

#ifndef BAYESIAN_NETWORKS_INFERENCE_REJECTION_SAMPLING_HPP
#define BAYESIAN_NETWORKS_INFERENCE_REJECTION_SAMPLING_HPP

#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <bayesian/network.hpp>
#include <bayesian/algorithm.hpp>
#include <bayesian/make_sample.hpp>

namespace bn {
namespace inference {

//! Provide a class performing Rejection Sampling for bn::network.
/*! @tparam RepresentMethod: Type of representation of network structure.
    TODO: Old-fashioned interface.
**/
template<class RepresentMethod>
class rejection_sampling {
public:
    using query_type = std::vector<component::random_variable_ptr>;
    using evidence_type = std::unordered_map<component::random_variable_ptr, std::size_t>;
    using sample_type = std::unordered_map<component::random_variable_ptr, matrix<std::size_t>>;
    using probability_type = std::unordered_map<component::random_variable_ptr, matrix<double>>;

    //! (ctor) Initialize the reasoner using the argument bayesian network.
    /*! @param[in]   network: bayesian network structure.
        @param[in]   cpt: CPT list corresponding to network.
    **/
    // TODO: tag
    rejection_sampling(network<RepresentMethod> const& network, cpt_manager const& cpts)
        : network_(network.clone()), cpts_(cpts)
    {
    }

    //! (ctor) Initialize the reasoner using the argument bayesian network and queries and evidences.
    /*! @param[in]   network: bayesian network structure.
        @param[in]   cpt: CPT list corresponding to network.
     **/
    // TODO: tag
    template<class NetworkType, class CptManagerType, class QueryType, class EvidenceType>
    rejection_sampling(NetworkType&& network, CptManagerType&& cpts, QueryType&& queries = QueryType(), EvidenceType&& evidences = EvidenceType())
        : network_(std::forward<NetworkType>(network)),
          cpts_(std::forward<CptManagerType>(cpts)),
          queries_(std::forward<QueryType>(queries)),
          evidences_(std::forward<EvidenceType>(evidences))
    {
    }
    
    //! Set a query variable, which is a target of reasoning.
    rejection_sampling& set_query(query_type queries)
    {
        if(!is_contained(queries, network_)) throw std::runtime_error(""); // TODO:
        queries_ = std::move(queries);
        is_modified_ = true;
        return *this;
    }

    //! Set query variables, which are targets of reasoning.
    template<class QueryType>
    rejection_sampling& set_query(QueryType&& queries)
    {
        if(!is_contained(queries, network_)) throw std::runtime_error(""); // TODO:
        queries_ = std::forward<QueryType>(queries);
        is_modified_ = true;
        return *this;
    }

    //! Set evidence variables, which are already known (conditions).
    rejection_sampling& set_evidence(evidence_type evidences)
    {
        evidences_ = std::move(evidences);
        is_modified_ = true;
        return *this;
    }

    //! Get query variables
    query_type const& query() const noexcept { return queries_; }

    //! Get query variables
    evidence_type const& evidence() const noexcept { return evidences_; }

    //! Get samples used in reasoning process.
    /*! Calling function `run` one or more times is required. */
    sample_type const& sample() const noexcept { return accepted_samples_; }

    //! Perform reasoning step in rejection sampling, which calculates query probabilities by using samples generated in sampling step.
    probability_type probability() const
    {
        probability_type probability_list;

        for(auto const& rv_sample : accepted_samples_)
        {
            auto const& raw_sample = rv_sample.second.data();
            auto const accepted_num = std::accumulate(raw_sample.cbegin(), raw_sample.cend(), 0);

            matrix<double> mat(rv_sample.second.sizes());
            for(std::size_t i = 0; i < raw_sample.size(); ++i)
                mat[std::vector<std::size_t>{i}] = static_cast<double>(raw_sample[i]) / accepted_num;

            probability_list.emplace(rv_sample.first, std::move(mat));
        }

        return probability_list;
    }

    // Re-initialize
    void reset()
    {
        initialize();
    }

    //! Perform sampling step in rejection sampling, until generating iterator_num
    sample_type run(std::size_t const iterator_num)
    {
        if(is_modified_)
        {
            initialize();
            is_modified_ = false;
        }

        auto working_accepted_samples = accepted_samples_;
        // std::vector<std::reference_wrapper<matrix<std::size_t>>> quick_accesser;
        // quick_accesser.reserve(queries_.size());
        // for(auto const& q : queries_)
        //     quick_accesser.push_back(working_accepted_samples[q]);

        for(std::size_t i = 0; i < iterator_num; ++i)
        {
            auto const sample = make_sample(network_, cpts_, topological_sorted_);
            if(!is_consistent(sample, evidences_)) continue;

            for(auto const& q : queries_)
                working_accepted_samples[q][std::vector<std::size_t>{sample.at(q)}] += 1;

            // for(std::size_t j = 0; j < quick_accesser.size(); ++j)
            //     quick_accesser[j][sample[queries_[j]]] += 1;
        }

        std::swap(working_accepted_samples, accepted_samples_);
        return accepted_samples_;
    }

private:
    // Refresh a generated sample and a topological order for graph.
    void initialize()
    {
        sample_type new_sample;
        for(auto const& q : queries_)
        {
            matrix<std::size_t> mat(std::vector<std::size_t>{ q->max_value }, 0);
            new_sample.emplace(q, std::move(mat));
        }

        auto topological_sorted = bn::topological_sort(network_);

        std::swap(new_sample, accepted_samples_);
        std::swap(topological_sorted, topological_sorted_);
    }

    // TODO: 外へ
    template<class Elem>
    bool is_contained(std::vector<Elem> const& elements, network<RepresentMethod> const& network) const
    {
        return std::all_of(elements.cbegin(), elements.cend(),
            [this, &network](auto const& element)
            {
                return is_contained(element, network);
            });
    }

    bool is_contained(component::random_variable_ptr const& rv, network<RepresentMethod> const& network) const
    {
        auto const& nodes = network_.all_node();
        return std::any_of(nodes.cbegin(), nodes.cend(),
            [this, &rv](auto const& n)
            {
                return n->get() == rv;
            });
    }

    bool is_contained(component::node_ptr const& node, network<RepresentMethod> const& network) const
    {
        auto const& nodes = network_.all_node();
        return std::any_of(nodes.cbegin(), nodes.cend(),
            [&node](auto const& n)
            {
                return n == node;
            });
    }

    // TODO: 外へ
    static bool is_consistent(cpt::condition_type const& sample, evidence_type const& evidence) // const
    {
        for(auto const& e : evidence)
        {
            if(sample.at(e.first) != e.second)
                return false;
        }

        return true;
    }

    bool is_modified_ = true;

    query_type queries_;
    evidence_type evidences_;

    sample_type accepted_samples_;

    network<RepresentMethod> network_;
    cpt_manager cpts_;
    std::vector<component::node_ptr> topological_sorted_;
};

} // namespace inference
} // namespace bn

#endif // BAYESIAN_NETWORKS_INFERENCE_REJECTION_SAMPLING_HPP

