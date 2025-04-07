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
 * Constructor for the HTwoHeuristic class.
 */
HTwoHeuristic::HTwoHeuristic(
    const shared_ptr<AbstractTask> &transform,
    bool cache_estimates, const string &description,
    utils::Verbosity verbosity)
    : Heuristic(transform, cache_estimates, description, verbosity),
      has_cond_effects(task_properties::has_conditional_effects(task_proxy)),
      goals(task_properties::get_fact_pairs(task_proxy.get_goals())) {
    if (log.is_at_least_normal()) {
        log << "Initializing h^2" << endl;
        log << "The implementation of the h^m heuristic is preliminary." << endl;
    }
    init_operator_caches();
}


/*
 * Computes the h^m value for a given state:
 */
int HTwoHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    }
    Tuple state_facts = task_properties::get_fact_pairs(state);
    init_hm_table(state_facts);
    init_operator_queue();
    update_hm_table();
    int h = eval(goals);
    if (h == INT_MAX) {
        return DEAD_END;
    }
    return h;
}

/**
* Sets up all auxiliary data structures concerning operators.
*/
void HTwoHeuristic::init_operator_caches() {
  	vector<int> empty_pre_op = {};
    for (auto op : task_proxy.get_operators()) {
        if (op.get_preconditions().empty()) {
            empty_pre_op.push_back(op.get_id());
        }
    }
    int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        int domain1_size = task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < domain1_size; ++j) {
			op_dict[FactPair(i, j)] = empty_pre_op;
        }
    }
    precondition_cache = {};
    partial_effect_cache = {};
    effect_conflict_cache.assign(task_proxy.get_operators().size(), std::vector<bool>(task_proxy.get_variables().size(), false));
	for (OperatorProxy op : task_proxy.get_operators()) {
        // Setup precondition cache
        Tuple preconditions = task_properties::get_fact_pairs(op.get_preconditions());
    	sort(preconditions.begin(), preconditions.end());
    	precondition_cache.push_back(preconditions);

        // Setup op_dict
        for (auto pre : preconditions) {
        	op_dict[pre].push_back(op.get_id());
        }

        // Check for operators without preconditions -> automatically add to op_dict
        if (preconditions.empty()) {
            empty_pre_op.push_back(op.get_id());
        }

        // Setup partial effect cache
    	Tuple effects;
    	for (EffectProxy eff : op.get_effects()) {
        	effects.push_back(eff.get_fact().get_pair());
            effect_conflict_cache[op.get_id()][eff.get_fact().get_pair().var] = true;
    	}
    	sort(effects.begin(), effects.end());
		partial_effect_cache.push_back(generate_all_pairs(effects));
    }
}

/*
 * Initializes h^m table.
 * If entry is contained in init state facts assigns 0, and infinity otherwise.
 * Pair containing variable -1 at second position indicates single fact.
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

/*
 * Check if Pair is contained in inital state facts. Unordered set allows constant time lookup.
 */
int HTwoHeuristic::check_in_initial_state(
    const Pair &hm_entry, const std::unordered_set<FactPair, FactPairHash> &state_facts_set) const {
    bool found_first = state_facts_set.find(hm_entry.first) != state_facts_set.end();
    bool found_second = (state_facts_set.find(hm_entry.second) != state_facts_set.end());
    return (found_first && found_second) ? 0 : INT_MAX;
}


void HTwoHeuristic::init_operator_queue() {
  	op_cost.assign(task_proxy.get_operators().size(), INT_MAX);
    changed_entries.assign(task_proxy.get_operators().size(), unordered_set<FactPair, FactPairHash>());
	for (OperatorProxy op : task_proxy.get_operators()) {
    	// Initialize operator queue with applicable operators
        if (is_op_applicable(precondition_cache[op.get_id()])) {
            op_queue.push_back(op.get_id());
            is_op_in_queue.insert(op.get_id());
        }
    }
}

/*
 * Check if op is applicable in initial state. Only works for initial state as it only considers single atom table entries.
 */
bool HTwoHeuristic::is_op_applicable(Tuple pre) const {
	for (auto fact : pre) {
    	if (hm_table.at(Pair(fact, FactPair(-1, -1))) != 0) {
        	return false;
        }
    }
    return true;
}


/*
 * Updates hm_table until no further improvements are made.
 */
void HTwoHeuristic::update_hm_table() {
    while (!op_queue.empty()) {
        OperatorProxy op = task_proxy.get_operators()[op_queue.front()];
        op_queue.pop_front();
        is_op_in_queue.erase(op.get_id());
        int c1 = eval(precondition_cache[op.get_id()]);
        if (c1 < op_cost[op.get_id()]) {
            changed_entries[op.get_id()].clear();
        	op_cost[op.get_id()] = c1;
        	for (Pair &partial_eff : partial_effect_cache[op.get_id()]) {
           		update_hm_entry(partial_eff, c1 + op.get_cost());
            	if (partial_eff.second.var == -1) {
                	extend_tuple(partial_eff.first, op, c1);
            	}
        	}
            continue;
        }
        if (c1 == INT_MAX) {
            continue;
        }
        handle_changed_entries(op);
    }
}


/*
 * Extends given partial effect by adding additional fact.
 */
void HTwoHeuristic::extend_tuple(const FactPair &f, const OperatorProxy &op, int eval) {
	Tuple pre = precondition_cache[op.get_id()];
    int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        if (effect_conflict_cache[op.get_id()][i]) {
        	continue;
        }
    	for (int j = 0; j < task_proxy.get_variables()[i].get_domain_size(); ++j) {
        	FactPair extend_fact = FactPair(i, j);
            // Check if extend_fact is reachable
            if (hm_table.at(Pair(extend_fact, FactPair(-1, -1))) == INT_MAX) {
            	continue;
            }
            Pair hm_pair = f.var > extend_fact.var ? Pair(extend_fact, f) : Pair(f, extend_fact);
            // Check if table entry can be updated with current op (without extend_Fact considered)
            if (hm_table.at(hm_pair) <= eval + op.get_cost()) {
            	continue;
            }
        	int c2 = extend_eval(extend_fact, pre, eval);
        	if (c2 != INT_MAX) {
            	update_hm_entry(hm_pair, c2 + op.get_cost());
        	}
        }
    }
}

/*
 * Handles changed entries for op. Used if op is achievable with same cost as in last iteration
 */
void HTwoHeuristic::handle_changed_entries(const OperatorProxy &op) {
	for (EffectProxy eff : op.get_effects()) {
    	FactPair effect = eff.get_fact().get_pair();
        for (FactPair entry : changed_entries[op.get_id()]) {
        	Pair hm_pair = effect.var > entry.var ? Pair(entry, effect) : Pair(effect, entry);
            if (hm_table.at(hm_pair) <= op_cost[op.get_id()] + op.get_cost()) {
            	continue;
            }
            int c2 = extend_eval(entry, precondition_cache[op.get_id()], op_cost[op.get_id()]);
        	if (c2 != INT_MAX) {
            	update_hm_entry(hm_pair, c2 + op.get_cost());
        	}
        }
    }
}

/*
 * Evaluates tuple by computing the maximum heuristic value among all its partial tuples. Used for pre(op) and goal.
 */
int HTwoHeuristic::eval(const Tuple &t) const {
    vector<Pair> pairs = generate_all_pairs(t);
    int max = 0;
    for (Pair &pair : pairs) {
        int h = hm_table.at(pair);
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
 * Evaluates extend_fact + pre. pre already evaluated with eval.
 */
int HTwoHeuristic::extend_eval(const FactPair &extend_fact, const Tuple &pre, int eval) const {
    int fact_eval = hm_table.at(Pair(extend_fact, FactPair(-1, -1)));
    int max = eval > fact_eval ? eval : fact_eval;
    for (FactPair fact0 : pre) {
      	if (fact0.var == extend_fact.var) {
            // Check if preconditions contradict extend_fact
          	if (fact0.value != extend_fact.value) {
                  return INT_MAX;
          	}
            // extend_fact ∈ pre
        	return eval;
        }
        Pair key = (fact0.var < extend_fact.var) ? Pair(fact0, extend_fact) : Pair(extend_fact, fact0);
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
 * Adds operators to the queue if f is precondition and was updated.
 */
void HTwoHeuristic::add_operator_to_queue(const Pair &p) {
    if (p.second.var == -1) {
    	for (int op_id : op_dict[p.first]) {
        	if (is_op_in_queue.find(op_id) == is_op_in_queue.end()) {
            	op_queue.push_back(op_id);
            	is_op_in_queue.insert(op_id);
        	}
    	}
        return;
    }

    for (int op_id : op_dict[p.first]) {
        if (is_op_in_queue.find(op_id) == is_op_in_queue.end()) {
            op_queue.push_back(op_id);
            is_op_in_queue.insert(op_id);
        }
        if (changed_entries[op_id].find(p.second) == changed_entries[op_id].end() && !effect_conflict_cache[op_id][p.second.var]) {
    		 changed_entries[op_id].insert(p.second);
    	}
    }
    for (int op_id : op_dict[p.second]) {
        if (is_op_in_queue.find(op_id) == is_op_in_queue.end()) {
            op_queue.push_back(op_id);
            is_op_in_queue.insert(op_id);
        }
        if (changed_entries[op_id].find(p.first) == changed_entries[op_id].end() && !effect_conflict_cache[op_id][p.first.var]) {
    		 changed_entries[op_id].insert(p.first);
    	}
    }

}

/*
 * Updates heuristic value of a pair in hm_table.
 * Affected operators are added to queue.
 */
int HTwoHeuristic::update_hm_entry(const Pair &p, int val) {
    if (hm_table[p] > val) {
        hm_table[p] = val;
        add_operator_to_queue(p);
        return val;
    }
    return -1;
}

/*
 * Generates all partial set of size <= 2 from given base tuple.
 */
vector<HTwoHeuristic::Pair> HTwoHeuristic::generate_all_pairs(const Tuple &base_tuple) const {
	vector<Pair> res;
    res.reserve(base_tuple.size() * (base_tuple.size() + 1) / 2);

    for (size_t i = 0; i < base_tuple.size(); ++i) {
        res.emplace_back(base_tuple[i], FactPair(-1, -1));

        for (size_t j = i + 1; j < base_tuple.size(); ++j) {
            res.emplace_back(base_tuple[i], base_tuple[j]);
        }
    }
    return res;
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
            get_heuristic_arguments_from_options(opts)
            );
    }
};

static plugins::FeaturePlugin<HTwoHeuristicFeature> _plugin;
}
