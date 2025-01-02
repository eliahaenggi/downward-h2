#ifndef HEURISTICS_HTWO_HEURISTIC_H
#define HEURISTICS_HTWO_HEURISTIC_H

#include "../heuristic.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <random>

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

    struct PairCompare {
        bool operator()(const std::pair<int, Pair>& p1, const std::pair<int, Pair>& p2) const {
            if (p1.first != p2.first) {
                return p1.first < p2.first;
            }
            if (p1.second.first != p2.second.first) {
                return p1.second.first < p2.second.first;
            }

            return p1.second.second < p2.second.second;
        }
    };

    struct CompareDistance {
        const std::unordered_map<Pair, int, PairHash>& distances;

        CompareDistance(const std::unordered_map<Pair, int, PairHash>& d) : distances(d) {}

        bool operator()(const Pair& a, const Pair& b) const {
            return distances.at(a) < distances.at(b);
        }
    };

    struct FactPairHash {
        size_t operator()(const FactPair &fact) const {
            return fact.var * 100003 + fact.value;
        }
    };

    struct OperatorInfo {

    	Tuple preconditions;
        std::vector<Pair> partial_effects;
        std::vector<FactPair> effects;

        OperatorInfo(const Tuple &pre, const std::vector<Pair> &par_eff, const std::vector<FactPair> eff) : preconditions(pre), partial_effects(par_eff), effects(eff) {}

    };

    // h^m table
    std::unordered_map<Pair, int, PairHash> hm_table;
    std::unordered_map<Pair, int, PairHash> distances;
    std::vector<Pair> table_order;
    std::unordered_map<FactPair, std::vector<OperatorProxy>, FactPairHash> operator_list;

    std::vector<OperatorInfo> operator_info_list;

    bool was_updated;

    // auxiliary methods
    void initialize_pairs(const Tuple &state_facts, std::vector<Pair> &init_pairs);
    void init_operator_info_list();
    void dijkstra(const std::vector<Pair> &init_pairs);

    void update_hm_table();
    int eval(const Tuple &t) const;
    int eval_preconditions(int &op_id) const;
    int hm_table_evaluation(const Tuple &t, const FactPair &fact, int eval) const;
    int update_hm_entry(const Pair &p, int val);
    void extend_tuple(const FactPair &f, const OperatorProxy &op, int eval);

    int check_in_initial_state(
        const Pair &hm_entry, const std::unordered_set<FactPair, FactPairHash> &state_facts_set) const;

    int get_operator_pre_value(const OperatorProxy &op, int var) const;
    bool contradict_effect_of(const OperatorProxy &op, int fact_var) const;

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
