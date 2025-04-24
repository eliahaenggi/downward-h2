#ifndef HEURISTICS_HTWO_HEURISTIC_H
#define HEURISTICS_HTWO_HEURISTIC_H

#include "../heuristic.h"

#include "../algorithms/priority_queues.h"

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

    // parameters
    const bool has_cond_effects;
    std::vector<FactPair> goals;


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
    int cost;
    int unsatisfied_preconditions;
    int id;

    BinaryOperator(std::vector<FactPair> pre,
                  Pair eff, int cost, int old_id, int i) : preconditions(pre), effect(eff), base_cost(cost), old_op_id(old_id), id(i) {}
	};
    // data structures
protected:
    std::unordered_map<Pair, int, PairHash> hm_table;
    std::unordered_map<Pair, std::vector<int>, PairHash> precondition_of;
    std::vector<BinaryOperator> binary_operators;
    priority_queues::AdaptiveQueue<Pair> queue;
    std::unordered_set<Pair, PairHash> partial_goals;

    void enqueue_if_necessary(Pair pair, int cost) {
        if (hm_table.at(pair) > cost) {
            hm_table[pair] = cost;
            queue.push(cost, pair);
        }
    }
    // Methods for initializing binary ops (called once)
    void init_binary_operators();
    void extend_binary_operators(const FactPair &f, const OperatorProxy &op, std::vector<bool>& effect_conflict);
    void setup_precondition_of();

    // Methods for initalizing hm_table (called once per eval)
    void init_hm_table(const std::vector<FactPair> &state_facts);
    int check_in_initial_state(
    const Pair &hm_entry, const std::unordered_set<FactPair, FactPairHash> &state_facts_set) const;

    // Methods for updating table
    int update_hm_table();

	std::vector<Pair> generate_all_pairs(const std::vector<FactPair> &base_tuple) const;

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
