#include "htwo_heuristic.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cassert>
#include <climits>
#include <set>

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
            << "It is rather slow." << endl
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
    init_operator_info_list();
    update_hm_table();
    int h = eval(goals);
    if (h == INT_MAX) {
        return DEAD_END;
    }
    return h;
}

/*
 * Initializes h^m table.
 * If tuple is contained in input tuple assigns 0, and infinity otherwise.
 */
void HTwoHeuristic::init_hm_table(const std::vector<FactPair> &state_facts) {
    unordered_set<FactPair, FactPairHash> state_facts_set(state_facts.begin(), state_facts.end());
    state_facts_set.insert(FactPair(-1, -1));

    int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        int domain1_size = task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < domain1_size; ++j) {
            Pair single_pair(FactPair(i, j), FactPair(-1, -1));
            hm_table[single_pair] = check_in_initial_state(single_pair, state_facts_set);
            for (int k = i + 1; k < num_variables; ++k) {
                int domain2_size = task_proxy.get_variables()[k].get_domain_size();
                for (int l = 0; l < domain2_size; ++l) {
                    Pair pair(FactPair(i, j), FactPair(k, l));
                    hm_table[pair] = check_in_initial_state(pair, state_facts_set);
                }
            }
        }
    }
}

int HTwoHeuristic::check_in_initial_state(
    const Pair &hm_entry, const std::unordered_set<FactPair, FactPairHash> &state_facts_set) const {
    bool found_first = state_facts_set.find(hm_entry.first) != state_facts_set.end();
    bool found_second = (state_facts_set.find(hm_entry.second) != state_facts_set.end());
    return (found_first && found_second) ? 0 : INT_MAX;
}

/**
* Generates partial effects (size <= 2) of all operators and saves them sorted in a map.
*/
void HTwoHeuristic::init_operator_info_list() {
	for (OperatorProxy op : task_proxy.get_operators()) {
        Tuple preconditions = task_properties::get_fact_pairs(op.get_preconditions());
    	sort(preconditions.begin(), preconditions.end());
    	Tuple effects;
    	for (EffectProxy eff : op.get_effects()) {
        	effects.push_back(eff.get_fact().get_pair());
    	}
    	sort(effects.begin(), effects.end());
		vector<Pair> partial_effs;
		generate_all_partial_tuples(effects, partial_effs);
        operator_info_list.push_back(OperatorInfo(preconditions, partial_effs));
    }
}


/*
 * Iteratively updates the h^m table until no further improvements are made.
 */
void HTwoHeuristic::update_hm_table() {
    do {
        was_updated = false;
        for (OperatorProxy op : task_proxy.get_operators()) {
            OperatorInfo &op_info = operator_info_list[op.get_id()];
            int c1 = eval(op_info.preconditions);
            if (c1 == INT_MAX) {
            	continue;
            }

            for (Pair &partial_eff : op_info.partial_effects) {
                update_hm_entry(partial_eff, c1 + op.get_cost());
                if (hm_table[partial_eff] != INT_MAX && partial_eff.second.var == -1) {

                    extend_tuple(partial_eff.first, op, c1);
                }
            }
        }
    } while (was_updated);
}


/*
 * Extends given partial effect by adding additional fact.
 */
void HTwoHeuristic::extend_tuple(const FactPair &f, const OperatorProxy &op, int eval) {
	Tuple pre = operator_info_list[op.get_id()].preconditions;
    int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        if (f.var == i || contradict_effect_of(op, i)) {
        	continue;
        }
    	for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
        	FactPair fact = FactPair(i, j);
        	int c2 = hm_table_evaluation(pre, fact, eval);
        	if (c2 != INT_MAX) {
            	Pair hm_pair = f.var > fact.var ? Pair(fact, f) : Pair(f, fact);
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
                  return INT_MAX;
            }
        	max = h;
        }
    }
    return max;
}

// Evaluates (t + fact). t-evaluation already given with eval.
int HTwoHeuristic::hm_table_evaluation(const Tuple &t, const FactPair &fact, int eval) const {
    int fact_eval = hm_table.at(Pair(fact, FactPair(-1, -1)));
    int max = eval > fact_eval ? eval : fact_eval;
    for (FactPair fact0 : t) {
      	if (fact0.var == fact.var) {
          	if (fact0.value != fact.value) {
                  return INT_MAX;
          	}
        	return eval;
        }
        Pair key = (fact0.var < fact.var) ? Pair(fact0, fact) : Pair(fact, fact0);
        int h = hm_table.at(key);

        if (h > max) {
        	if (h == INT_MAX) {
                return INT_MAX;
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

bool HTwoHeuristic::contradict_effect_of(
    const OperatorProxy &op, int fact_var) const {
    for (EffectProxy eff : op.get_effects()) {
        if (eff.get_fact().get_variable().get_id() == fact_var) {
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
        res.emplace_back(Pair(base_tuple[i], FactPair(-1, -1)));

        for (size_t j = i + 1; j < base_tuple.size(); ++j) {
            res.emplace_back(Pair(base_tuple[i], base_tuple[j]));
        }
    }
}

void HTwoHeuristic::print_table() const {
    stringstream ss;
    for (auto entry : hm_table) {
      	if (entry.second == INT_MAX) {
        	continue;
      	}
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
