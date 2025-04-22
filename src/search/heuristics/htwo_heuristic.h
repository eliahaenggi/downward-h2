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
class HTwoHeuristic : public Heuristic {
    protected:
    using Tuple = std::vector<FactPair>;

    // parameters
    const bool has_cond_effects;
    Tuple goals;


    struct Pair {
        FactPair first;
        FactPair second;
        std::size_t hash;

        Pair(const FactPair &f, const FactPair &s) : first(f), second(s), hash(compute_hash(first, second)) {}

        bool operator==(const Pair &other) const {
            return first == other.first && second == other.second;
        }

    protected:
        static std::size_t compute_hash(const FactPair &f1, const FactPair &f2) {
            const int MOD = 100003; // Prime
            std::size_t h1 = f1.var * MOD + f1.value;
            std::size_t h2 = f2.var * MOD + f2.value;
            return h2 * MOD + h1;
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


    struct BinaryOperator {
    std::vector<FactPair> preconditions;
    Pair effect;
    int base_cost;
    int old_op_id;

    BinaryOperator(std::vector<FactPair> pre,
                  Pair eff, int cost, int id) : preconditions(pre), effect(eff), base_cost(cost), old_op_id(id) {}
	};
    // data structures
protected:
    std::unordered_map<Pair, int, PairHash> hm_table;
    std::vector<BinaryOperator> binary_operators;
    std::deque<int> op_queue;
    bool was_updated;

    // Auxiliary data structurs that speed up implementation (Could also be removed in case of memory issues)
    std::unordered_set<int> is_op_in_queue; // stores all operators that are in queue for constant time look up
    std::vector<Tuple> precondition_cache;
    std::vector<std::vector<Pair>> partial_effect_cache;
    std::vector<std::vector<bool>> effect_conflict_cache; // Stores if variable is in effect of operator

   	std::vector<int> op_cost;
    std::vector<std::unordered_set<FactPair, FactPairHash>> changed_entries;
    // Stores for each FactPair a list of operators where the fact occures in pre
    mutable std::unordered_map<FactPair, std::vector<int>, FactPairHash> op_dict;

    // Methods for initalizing data structures
    void init_hm_table(const Tuple &state_facts);
    int check_in_initial_state(
    const Pair &hm_entry, const std::unordered_set<FactPair, FactPairHash> &state_facts_set) const;
    void init_operator_caches();
    void extend_binary_operators(const FactPair &f, const OperatorProxy &op);
    void init_operator_queue();
    bool is_op_applicable(Tuple pre) const;

    // Methods for updating table
    void update_hm_table();
    void extend_tuple(const FactPair &f, const OperatorProxy &op, int eval);
    int eval(const Tuple &t) const;
    inline int extend_eval(const FactPair &extend_fact, const Tuple &pre, int eval) const;
    inline void extend_changed_entry(const OperatorProxy &op);

    inline void update_hm_entry(const Pair &p, int val);
    inline void add_operator_to_queue(const Pair &p);

	std::vector<Pair> generate_all_pairs(const Tuple &base_tuple) const;

    void print_table() const;

protected:
    virtual int compute_heuristic(const State &ancestor_state) override;

public:
    HTwoHeuristic(
        const std::shared_ptr<AbstractTask> &transform,
        bool cache_estimates, const std::string &description,
        utils::Verbosity verbosity);

};
}

#endif
