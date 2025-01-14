#ifndef HEURISTICS_HTWO_HEURISTIC_H
#define HEURISTICS_HTWO_HEURISTIC_H

#include "../heuristic.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <deque>



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

    // parameters
    const int m;
    const bool has_cond_effects;

    const Tuple goals;


    struct Pair {
        FactPair first;
        FactPair second;
        std::size_t hash;

        Pair(const FactPair &f, const FactPair &s) : first(f), second(s), hash(compute_hash(first, second)) {}

        bool operator==(const Pair &other) const {
            return first == other.first && second == other.second;
        }

    private:
        static std::size_t compute_hash(const FactPair &f1, const FactPair &f2) {
            const int MOD = 100003; // Prime
            std::size_t h1 = f1.var * MOD + f1.value;
            std::size_t h2 = f2.var * MOD + f2.value;
            return h1 * MOD + h2;
        }
    };

    struct PairHash {
        std::size_t operator()(const Pair &pair) const {
            return pair.hash;
        }
    };

    struct FactPairHash {
        size_t operator()(const FactPair &fact) const {
            return fact.var * 100003 + fact.value;
        }
    };

    // data structures
    std::unordered_map<Pair, int, PairHash> hm_table;

    std::vector<Tuple> precondition_cache;
    std::vector<std::vector<Pair>> partial_effect_cache;
    mutable std::unordered_map<FactPair, std::vector<OperatorProxy>, FactPairHash> op_dictionary;

	std::deque<OperatorProxy> op_queue;
    std::unordered_set<int> is_op_in_queue; // for constant check if op in queue

    // Methods for initalizing data structures
    void init_hm_table(const Tuple &state_facts);
    int check_in_initial_state(
    const Pair &hm_entry, const std::unordered_set<FactPair, FactPairHash> &state_facts_set) const;
    void init_operator_caches();
    bool is_op_applicable(Tuple pre) const;

    // Methods for updating table
    void update_hm_table();
    void extend_tuple(const FactPair &f, const OperatorProxy &op, int eval);
    int eval(const Tuple &t) const;
    int extend_eval(const FactPair &extend_fact, const Tuple &pre, int eval) const;

    int update_hm_entry(const Pair &p, int val);
    void add_operator_to_queue(const FactPair &f);

    bool contradict_effect_of(const OperatorProxy &op, int fact_var) const;

	std::vector<Pair> generate_all_pairs(const Tuple &base_tuple) const;

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
