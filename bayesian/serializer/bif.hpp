#ifndef BNI_SERIALIZER_BIF_HPP
#define BNI_SERIALIZER_BIF_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <boost/spirit/home/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <bayesian/graph.hpp>

namespace bn {
namespace serializer {

// 辛い
namespace qi  = boost::spirit::qi;
namespace phx = boost::phoenix;

class bif {
public:
    // パース直後の木を形成する
    struct variable_holder {
        std::string variable_name;
        unsigned int possible_value;
        std::vector<std::string> possible_value_name;
    };

    struct network_holder {
        std::string network_name;
    };

    struct probability_holder {
        typedef std::pair<std::string, std::vector<std::string>> rv_condition_type;
        typedef std::vector<std::pair<std::vector<std::string>, std::vector<double>>> probability_type;
        rv_condition_type rv_condition;
        probability_type probability;
    };

    template<class Iterator>
    bn::graph_t parse(Iterator const& begin, Iterator const& end)
    {
        parse_grammar<Iterator> grammar;

        // Parse
        auto it = begin;
        bool success = qi::phrase_parse(it, end, grammar, qi::ascii::space);
        if(!success || it != end)
        {
            // 木構造作成が正常でない(パース失敗)
            throw std::runtime_error("cannot parse bif");
        }

        //
        // 木構造をグラフとして作成していく
        //
        bn::graph_t graph;
        std::unordered_map<std::string, std::pair<bn::vertex_type, variable_holder>> dictionary;

        // ノードについての作成・辞書作成
        for(auto const& node : grammar.variables)
        {
            auto vertex = graph.add_vertex();
            vertex->selectable_num = node.possible_value;
            dictionary[node.variable_name] = std::pair<bn::vertex_type, variable_holder>(vertex, node);
        }

        // CPTの読み取り
        for(auto const& probability : grammar.probabilities)
        {
            // 注目しているノードについて
            auto const& random_variables = probability.rv_condition;
            auto const& target_vertex = dictionary[random_variables.first].first;

            // ノードの親についてまとめ
            // CPTの雛形まで作成する
            std::vector<bn::vertex_type> conditional_vertexes;
            for(auto const& conditional_node : random_variables.second)
            {
                auto const& conditional_vertex = dictionary[conditional_node].first;
                graph.add_edge(conditional_vertex, target_vertex);  // 辺を張る
                conditional_vertexes.push_back(conditional_vertex); // target_vertexの条件として保存
            }
            target_vertex->cpt.assign(conditional_vertexes, target_vertex); // 雛形

            // CPTが指定されていれば，CPTを確定
            // されていなければループに入らない
            for(auto const& one_condition : probability.probability)
            {
                // 条件として与えられた文字列を便宜的に数字(index)に直し，条件とする．
                bn::condition_t cond;
                std::transform(
                    one_condition.first.begin(), one_condition.first.end(),
                    random_variables.second.begin(),
                    std::inserter(cond, cond.begin()),
                    [&dictionary](std::string const& value, std::string const& label) -> bn::condition_t::value_type
                    {
                        auto const& cond_node = dictionary[label];
                        auto const index = std::distance(
                            cond_node.second.possible_value_name.begin(),
                            std::find(cond_node.second.possible_value_name.begin(), cond_node.second.possible_value_name.end(), value)
                            );

                        return std::make_pair(cond_node.first, static_cast<int>(index));
                    });

                target_vertex->cpt[cond].second = one_condition.second;
            }
        }

        return graph;
    }

private:
    template<class Iterator>
    struct parse_grammar : qi::grammar<Iterator, void(), qi::ascii::space_type> {
        //
        // Functions
        //
        std::function<std::pair<std::vector<std::string>, std::vector<double>>(std::vector<double>)> nonconditional;

        //
        // Rule of "variable" Section
        //
        qi::rule<Iterator, std::vector<std::string>(), qi::ascii::space_type> enumeration_rule;
        qi::rule<Iterator, unsigned int(), qi::ascii::space_type> array_size_rule;
        qi::rule<Iterator, variable_holder(), qi::ascii::space_type> variable_rule;

        //
        // Rule of "network" Section
        //
        qi::rule<Iterator, network_holder(), qi::ascii::space_type> network_rule;

        //
        // Rule of "probability" Section
        //
        qi::rule<Iterator, std::pair<std::string, std::vector<std::string>>(), qi::ascii::space_type> rv_enum_rule;
        qi::rule<Iterator, std::vector<std::string>(), qi::ascii::space_type> condition_enum_rule;
        qi::rule<Iterator, std::vector<double>(), qi::ascii::space_type> probability_enum_rule;
        qi::rule<Iterator, std::pair<std::vector<std::string>, std::vector<double>>(), qi::ascii::space_type> probability_line_rule;
        qi::rule<Iterator, std::pair<std::vector<std::string>, std::vector<double>>(), qi::ascii::space_type> independent_line_rule;
        qi::rule<Iterator, probability_holder(), qi::ascii::space_type> rv_rule;
        qi::rule<Iterator, probability_holder(), qi::ascii::space_type> probability_rule;

        //
        // Rule of ALL
        //
        qi::rule<Iterator, void(), qi::ascii::space_type> global_rule;

        //
        // Parsed Holder
        //
        std::vector<network_holder> networks;
        std::vector<variable_holder> variables;
        std::vector<probability_holder> probabilities;

        parse_grammar() : parse_grammar::base_type(global_rule)
        {
            //
            // Functions
            //
            nonconditional =
                [](std::vector<double> target)
                {
                    return std::make_pair(std::vector<std::string>(), std::move(target));
                };

            //
            // Rule of "variable" Section
            //
            enumeration_rule =
                qi::lit("{")
                >> +qi::alnum % qi::lit(",")
                >> qi::lit("}");
            array_size_rule =
                qi::lit("[") >> qi::uint_ >> qi::lit("]");
            variable_rule =
                qi::lit("variable")
                >> +qi::alnum
                >> qi::lit("{")
                >> qi::lit("type") >> qi::lit("discrete")>> array_size_rule >> enumeration_rule >> qi::lit(";")
                >> qi::lit("}");

            //
            // Rule of "network" Section
            //
            network_rule =
                qi::lit("network")
                >> +qi::alnum
                >> qi::lit("{") >> qi::lit("}");

            //
            // Rule of "probability" Section
            //
            rv_enum_rule =
                qi::lit("(")
                >> (+qi::alnum)
                >>-( +qi::lit("|")
                >> (+qi::alnum) % qi::lit(",") )
                >> qi::lit(")");
            condition_enum_rule =
                qi::lit("(")
                >> +qi::alnum % qi::lit(",")
                >> qi::lit(")");
            probability_enum_rule =
                qi::double_ % qi::lit(",");
            probability_line_rule =
                condition_enum_rule >> probability_enum_rule >> qi::lit(";");
            independent_line_rule =
                qi::lit("table") >> probability_enum_rule[qi::_val = phx::bind(nonconditional, qi::_1)] >> qi::lit(";");
            rv_rule =
                rv_enum_rule >> qi::lit("{")
                >> -(+independent_line_rule | +probability_line_rule)
                >> qi::lit("}");
            probability_rule =
                qi::lit("probability") >> rv_rule;

            //
            // Rule of ALL
            //
            global_rule =
                *(
                    network_rule[phx::push_back(phx::ref(networks), qi::_1)]
                    | variable_rule[phx::push_back(phx::ref(variables), qi::_1)]
                    | probability_rule[phx::push_back(phx::ref(probabilities), qi::_1)]
                );

        }
    };

};

} // namespace serializer
} // namespace bn

BOOST_FUSION_ADAPT_STRUCT(
    bn::serializer::bif::variable_holder,
    (std::string, variable_name)
    (unsigned int, possible_value)
    (std::vector<std::string>, possible_value_name)
)

BOOST_FUSION_ADAPT_STRUCT(
    bn::serializer::bif::network_holder,
    (std::string, network_name)
)

BOOST_FUSION_ADAPT_STRUCT(
    bn::serializer::bif::probability_holder,
    (bn::serializer::bif::probability_holder::rv_condition_type, rv_condition)
    (bn::serializer::bif::probability_holder::probability_type, probability)
)

#endif // #ifndef BNI_SERIALIZER_BIF_HPP

