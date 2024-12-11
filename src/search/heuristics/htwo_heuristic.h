#ifndef HEURISTICS_HTWO_HEURISTIC_H
#define HEURISTICS_HTWO_HEURISTIC_H

#include "../heuristic.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

namespace plugins {
class Options;
}

namespace htwo_heuristic {
/*
  Haslum's h^m heuristic family ("critical path heuristics").

  This is a very slow implementation and should not be used for
  speed benchmarks.
*/

class HTwoHeuristic : public Heuristic {
    using Tuple = std::vector<FactPair>;
    using Pair = std::pair<FactPair, FactPair>;

    // parameters
    const int m;
    const bool has_cond_effects;

    const Tuple goals;


    struct PairHash {
        std::size_t operator()(const Pair& pair) const {
            std::string hash = std::to_string(pair.first.var) + std::to_string(pair.first.value);
            hash += std::to_string(pair.second.var) + std::to_string(pair.second.value);
            return std::hash<std::string>()(hash);
        }
    };
    // h^m table
    std::unordered_map<Pair, int, PairHash> hm_table;


    mutable std::unordered_map<int, Tuple> precondition_cache;

    bool was_updated;

    // auxiliary methods
    void init_hm_table(const Tuple &state_facts);
    void update_hm_table();
    int eval(const Tuple &t) const;
    int update_hm_entry(const Pair &p, int val);
    void extend_tuple(const Pair &p, const OperatorProxy &op);

    bool check_tuple_in_tuple(const Pair &tuple, const Tuple &big_tuple) const;

    int get_operator_pre_value(const OperatorProxy &op, int var) const;
    Tuple get_operator_pre(const OperatorProxy &op) const;
    Tuple get_operator_eff(const OperatorProxy &op) const;
    bool contradict_effect_of(const OperatorProxy &op, FactPair fact) const;

	void generate_all_partial_tuples(const Tuple &base_tuple,
                                     std::vector<Pair> &res) const;

    void print_table() const;

protected:
    virtual int compute_heuristic(const State &ancestor_state) override;

public:
    HTwoHeuristic(
        int m, const std::shared_ptr<AbstractTask> &transform,
        bool cache_estimates, const std::string &description,
        utils::Verbosity verbosity);

    virtual bool dead_ends_are_reliable() const override;
};
}

#endif
