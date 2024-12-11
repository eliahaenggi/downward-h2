#include "htwo_heuristic.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cassert>
#include <limits>
#include <set>
#include <thread>
#include <chrono>

using namespace std;

namespace htwo_heuristic {
/*
 * Constructor for the HMHeuristic class.
 * Precomputes all possible tuples of size <= m.
 */
HTwoHeuristic::HTwoHeuristic(
    int m, const shared_ptr<AbstractTask> &transform,
    bool cache_estimates, const string &description,
    utils::Verbosity verbosity)
    : Heuristic(transform, cache_estimates, description, verbosity),
      m(m),
      has_cond_effects(task_properties::has_conditional_effects(task_proxy)),
      goals(task_properties::get_fact_pairs(task_proxy.get_goals())) {
    if (log.is_at_least_normal()) {
        log << "Using h^" << m << "." << endl;
        log << "The implementation of the h^m heuristic is preliminary." << endl
            << "It is SLOOOOOOOOOOOW." << endl
            << "Please do not use this for comparison!" << endl;
    }
}


bool HTwoHeuristic::dead_ends_are_reliable() const {
    return !task_properties::has_axioms(task_proxy) && !has_cond_effects;
}


/*
 * Computes the h^m value for a given state:
 * Checks if state is a goal state (heuristic = 0 if true). Initializes h^m table with state facts.
 * Updates h^m table to propagate values. Evaluates goal facts to compute the heuristic value.
 */
int HTwoHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    }
    Tuple state_facts = task_properties::get_fact_pairs(state);

    init_hm_table(state_facts);
    update_hm_table();
    int h = eval(goals);

    if (h == numeric_limits<int>::max()) {
        return DEAD_END;
    }
    return h;
}

/*
 * Initializes h^m table.
 * If tuple is contained in input tuple assigns 0, and infinity otherwise.
 */
void HTwoHeuristic::init_hm_table(const Tuple &state_facts) {
	int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
      	int domain1_size = task_proxy.get_variables()[i].get_domain_size();
    	for (int j = 0; j < domain1_size; ++j) {
            Pair single_pair(FactPair(i, j), FactPair(-1, -1));
            hm_table[single_pair] = check_tuple_in_tuple(single_pair, state_facts) ? 0 : numeric_limits<int>::max();
            for (int k = i + 1; k < num_variables; ++k) {
            	int domain2_size = task_proxy.get_variables()[k].get_domain_size();
            	for (int l = 0; l < domain2_size; ++l) {
                  Pair pair(FactPair(i, j), FactPair(k, l));
                  hm_table[pair] = check_tuple_in_tuple(single_pair, state_facts) ? 0 : numeric_limits<int>::max();
            	}
            }
    	}
    }
}

/*
 * Iteratively updates the h^m table until no further improvements are made.
 * For each operator:
 * - Computes the h value of preconditions.
 * - If preconditions are reachable generates all partial effect tuples for table update.
 * - Extends partial tuples with size < m to explore potential improvements.
 */
void HTwoHeuristic::update_hm_table() {
    do {
        was_updated = false;

        for (OperatorProxy op : task_proxy.get_operators()) {
            Tuple pre = get_operator_pre(op);

            int c1 = eval(pre);
            if (c1 != numeric_limits<int>::max()) {
                Tuple eff = get_operator_eff(op);
                vector<Pair> partial_effs;
                generate_all_partial_tuples(eff, partial_effs);
                for (Pair &partial_eff : partial_effs) {
                    update_hm_entry(partial_eff, c1 + op.get_cost());

                    if (partial_eff.second.var == -1 ) {
                        extend_tuple(partial_eff, op);
                    }
                }
            }
        }
    } while (was_updated);
}


/*
 * Extends given tuple by adding additional facts.
 * Checks for contradictions between operator effects and tuple.
 * If no contradiction exists, updates h^m table if improvements are found.
 */
void HTwoHeuristic::extend_tuple(const Pair &p, const OperatorProxy &op) {
    for (const auto &hm_ent : hm_table) {
        const Pair &hm_pair = hm_ent.first;

        if (hm_pair.second.var == -1) {
            continue;
        }

        if (contradict_effect_of(op, hm_pair.first) || contradict_effect_of(op, hm_pair.second)) {
            continue;
        }

        FactPair fact = FactPair(-1, -1);
        if (p.first == hm_pair.first) {
            fact = hm_pair.second;
        } else if (p.first == hm_pair.second) {
            fact = hm_pair.first;
        }

        if (fact.var == -1) {
            continue;
        }

        Tuple pre = get_operator_pre(op);
        auto it = std::lower_bound(pre.begin(), pre.end(), fact);
        if (it == pre.end() || *it != fact) {
            pre.insert(it, fact);
        }

        std::unordered_set<int> vars;
        bool is_valid = true;
        for (const FactPair &f : pre) {
            if (!vars.insert(f.var).second) {
                is_valid = false;
                break;
            }
        }

        // Update the heuristic table if valid
        if (is_valid) {
            int c2 = eval(pre);
            if (c2 != std::numeric_limits<int>::max()) {
                update_hm_entry(hm_pair, c2 + op.get_cost());
            }
        }
    }
}



/*
 * Evaluates tuple by computing the maximum heuristic value among all its partial tuples.
 */
int HTwoHeuristic::eval(const Tuple &t) const {
    vector<Pair> partial;
    generate_all_partial_tuples(t, partial);
    int max = 0;
    for (Pair &pair : partial) {
        int h = hm_table.at(pair);

        if (h > max) {
        	if (h == numeric_limits<int>::max()) {
                  return numeric_limits<int>::max();
            }
        	max = h;
        }

    }
    return max;
}

/*
 * Updates the heuristic value of a tuple in the h^m table.
 * Sets "was_updated" flag to true to indicate a change occurred.
 */
int HTwoHeuristic::update_hm_entry(const Pair &p, int val) {
    if (hm_table[p] > val) {
        hm_table[p] = val;
        was_updated = true;
    }
    return val;
}

bool HTwoHeuristic::check_tuple_in_tuple(
    const Pair &pair, const Tuple &big_tuple) const {
    bool found_first = false;
    bool found_second = pair.second.var == -1;
    for (auto &fact : big_tuple) {
        if (!found_first && fact == pair.first) {
            found_first = true;
        } else if (!found_second && fact == pair.second) {
            found_second = true;
        }
        if (found_first && found_second) {
            return true;
        }
    }
    return false;
}


HTwoHeuristic::Tuple HTwoHeuristic::get_operator_pre(const OperatorProxy &op) const {
    int op_id = op.get_id();

    auto it = precondition_cache.find(op_id);
    if (it != precondition_cache.end()) {
        return it->second;
    }

    Tuple preconditions = task_properties::get_fact_pairs(op.get_preconditions());
    std::sort(preconditions.begin(), preconditions.end());
    precondition_cache[op_id] = preconditions;

    return preconditions;
}


HTwoHeuristic::Tuple HTwoHeuristic::get_operator_eff(const OperatorProxy &op) const {
    Tuple effects;
    for (EffectProxy eff : op.get_effects()) {
        effects.push_back(eff.get_fact().get_pair());
    }
    return effects;
}


bool HTwoHeuristic::contradict_effect_of(
    const OperatorProxy &op, FactPair fact) const {
    for (EffectProxy eff : op.get_effects()) {
        FactProxy eff_fact = eff.get_fact();
        if (eff_fact.get_variable().get_id() == fact.var && eff_fact.get_value() != fact.value) {
            return true;
        }
    }
    return false;
}

/*
 * Generates all partial tuples of size <= m from given base tuple.
 */
void HTwoHeuristic::generate_all_partial_tuples(
    const Tuple &base_tuple, vector<Pair> &res) const {
    res.reserve(base_tuple.size() * (base_tuple.size() + 1) / 2);

    for (size_t i = 0; i < base_tuple.size(); ++i) {
        res.emplace_back(base_tuple[i], FactPair(-1, -1));

        for (size_t j = i + 1; j < base_tuple.size(); ++j) {
            res.emplace_back(base_tuple[i], base_tuple[j]);
        }
    }
}

void HTwoHeuristic::print_table() const {
    stringstream ss;
    for (auto entry : hm_table) {
        Pair pair = entry.first;
        ss << "[" << pair.first.var << " = " << pair.first.value << ", " << pair.second.var << " = " << pair.second.value << "] = " << entry.second << endl;
    }
    log << ss.str() << endl;
}

class HTwoHeuristicFeature
    : public plugins::TypedFeature<Evaluator, HTwoHeuristic> {
public:
    HTwoHeuristicFeature() : TypedFeature("h2") {
        document_title("h^2 heuristic");

        add_option<int>("m", "subset size", "2", plugins::Bounds("1", "infinity"));
        add_heuristic_options_to_feature(*this, "h2");

        document_language_support("action costs", "supported");
        document_language_support("conditional effects", "ignored");
        document_language_support("axioms", "ignored");

        document_property(
            "admissible",
            "yes for tasks without conditional effects or axioms");
        document_property(
            "consistent",
            "yes for tasks without conditional effects or axioms");
        document_property(
            "safe",
            "yes for tasks without conditional effects or axioms");
        document_property("preferred operators", "no");
    }

    virtual shared_ptr<HTwoHeuristic> create_component(
        const plugins::Options &opts,
        const utils::Context &) const override {
        return plugins::make_shared_from_arg_tuples<HTwoHeuristic>(
            opts.get<int>("m"),
            get_heuristic_arguments_from_options(opts)
            );
    }
};

static plugins::FeaturePlugin<HTwoHeuristicFeature> _plugin;
}
