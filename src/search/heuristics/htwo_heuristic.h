#ifndef HEURISTICS_HTWO_HEURISTIC_H
#define HEURISTICS_HTWO_HEURISTIC_H

#include "../heuristic.h"

#include <algorithm>
#include <iostream>
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
    std::size_t operator()(const Pair &pair) const {
        const int MOD = 1009; // Prime greater than max domain size
        std::size_t h1 = pair.first.var * MOD + pair.first.value;
        std::size_t h2 = pair.second.var * MOD + pair.second.value;
        return h1 * MOD + h2;
    }
};
    // h^m table
    std::unordered_map<Pair, int, PairHash> hm_table;
    mutable std::unordered_map<int, Tuple> precondition_cache;
    mutable std::unordered_map<int, std::vector<Pair>> partial_effect_cache;

    bool was_updated;

    // auxiliary methods
    void init_hm_table(const Tuple &state_facts);
    void init_partial_effects();
    void update_hm_table();
    int eval(const Tuple &t) const;
    int hm_table_evaluation(const Tuple &t, const FactPair &fact, int eval) const;
    int update_hm_entry(const Pair &p, int val);
    void extend_tuple(const Pair &p, const OperatorProxy &op, int eval);

    bool check_in_initial_state(const Pair &hm_entry, const Tuple &state_facts) const;

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
