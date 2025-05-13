#include "htwo_heuristic.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"

#include <climits>

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
    vector<FactPair> state_facts = task_properties::get_fact_pairs(state);
    init_hm_table(state_facts);
    init_operator_queue();
    update_hm_table();
    int h = eval(goals);
    if (h == INT_MAX) {
        return DEAD_END;
    }
    return h;
}

/*
* Sets up all auxiliary operator data structures. Since these do not change throughout the search process, this is done once in the constructor.
*/
void HTwoHeuristic::init_operator_caches() {
  	vector<int> empty_pre_op = {};
    for (auto op : task_proxy.get_operators()) {
        if (op.get_preconditions().empty()) {
            empty_pre_op.push_back(op.get_id());
        }
    }
    const int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        const int domain1_size = task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < domain1_size; ++j) {
			op_dict[FactPair(i, j)] = empty_pre_op;
        }
    }
    precondition_cache = {};
    partial_effect_cache = {};
    effect_conflict_cache.assign(task_proxy.get_operators().size(), std::vector<bool>(task_proxy.get_variables().size(), false));
	for (OperatorProxy op : task_proxy.get_operators()) {
        // Setup precondition cache
        vector<FactPair> preconditions = task_properties::get_fact_pairs(op.get_preconditions());
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
        vector<FactPair> effects;
    	for (EffectProxy eff : op.get_effects()) {
        	effects.push_back(eff.get_fact().get_pair());
            effect_conflict_cache[op.get_id()][eff.get_fact().get_pair().var] = true;
    	}
    	sort(effects.begin(), effects.end());
		partial_effect_cache.push_back(generate_all_pairs(effects));
    }
}

/*
 * Initializes hm table.
 * If entry is contained in init_state_atoms assigns 0, and infinity otherwise.
 * Pair containing variable -1 at second position indicates single fact.
 */
void HTwoHeuristic::init_hm_table(const std::vector<FactPair> &init_state_atoms) {
    unordered_set<FactPair, FactPairHash> state_atoms_set(init_state_atoms.begin(), init_state_atoms.end());
        state_atoms_set.insert(FactPair(-1, -1));

    const int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        const int domain1_size = task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < domain1_size; ++j) {
            Pair single_pair(FactPair(i, j), FactPair(-1, -1));
            hm_table[single_pair] = check_in_initial_state(single_pair, state_atoms_set);

            for (int k = i + 1; k < num_variables; ++k) {
                const int domain2_size = task_proxy.get_variables()[k].get_domain_size();
                for (int l = 0; l < domain2_size; ++l) {
                    Pair pair(FactPair(i, j), FactPair(k, l));
                    hm_table[pair] = check_in_initial_state(pair, state_atoms_set);
                }
            }
        }
    }
}

/*
 * Check if Pair is contained in initial state atoms. Unordered set allows constant time lookup.
 */
int HTwoHeuristic::check_in_initial_state(
    const Pair &hm_entry, const std::unordered_set<FactPair, FactPairHash> &state_atoms_set) const {
    bool found_first = state_atoms_set.find(hm_entry.first) != state_atoms_set.end();
    bool found_second = (state_atoms_set.find(hm_entry.second) != state_atoms_set.end());
    return (found_first && found_second) ? 0 : INT_MAX;
}

/*
 * Adds operators to queue if they are applicable in the initial state.
 */
void HTwoHeuristic::init_operator_queue() {
  	op_cost.assign(task_proxy.get_operators().size(), INT_MAX);
    changed_entries.assign(task_proxy.get_operators().size(), unordered_set<FactPair, FactPairHash>());
	for (OperatorProxy op : task_proxy.get_operators()) {
        const int op_id = op.get_id();
    	// Initialize operator queue with applicable operators
        if (is_op_applicable(precondition_cache[op_id])) {
            op_queue.push_back(op_id);
            is_op_in_queue.insert(op_id);
        }
    }
}

/*
 * Check if op is applicable in initial state. Only works for initial state as it only considers single atom table entries.
 * Precondition-free operators are always applicable.
 */
bool HTwoHeuristic::is_op_applicable(vector<FactPair> pre) const {
	for (auto fact : pre) {
    	if (hm_table.find(Pair(fact, FactPair(-1, -1)))->second != 0) {
        	return false;
        }
    }
    return true;
}


/*
 * Updates hm_table until operator queue is empty and with it, no further improvements are possible.
 */
void HTwoHeuristic::update_hm_table() {
    while (!op_queue.empty()) {
        const int op_id = op_queue.front();
        OperatorProxy op = task_proxy.get_operators()[op_id];
        const int cost = op.get_cost();
        op_queue.pop_front();
        is_op_in_queue.erase(op_id);
        int c1 = eval(precondition_cache[op_id]);
        if (c1 == op_cost[op_id]) {
          	if (c1 != INT_MAX) {
		    	extend_changed_entry(op);
        	}
            continue;
        }
        changed_entries[op_id].clear();
        op_cost[op_id] = c1;
        for (Pair &partial_eff : partial_effect_cache[op_id]) {
        	update_hm_entry(partial_eff, c1 + cost);
            if (partial_eff.second.var == -1) {
                extend_entry(partial_eff.first, op, c1);
            }
        }
    }
}


/*
 * Extends given partial effect by adding additional atom.
 */
void HTwoHeuristic::extend_entry(const FactPair &f, const OperatorProxy &op, int eval) {
    const auto &variables = task_proxy.get_variables();
    const int op_id = op.get_id();
    const int op_cost = op.get_cost();
    const vector<FactPair> &pre = precondition_cache[op_id];
    for (size_t i = 0; i < variables.size(); ++i) {
        if (effect_conflict_cache[op_id][i]) {
        	continue;
        }
        const int domain_size = variables[i].get_domain_size();
    	for (int j = 0; j < domain_size; ++j) {
        	const FactPair extend_fact = FactPair(i, j);
            // Check if extend_fact is reachable
            if (hm_table.find(Pair(extend_fact, FactPair(-1, -1)))->second == INT_MAX) {
            	continue;
            }
            const Pair hm_pair = f.var > extend_fact.var ? Pair(extend_fact, f) : Pair(f, extend_fact);
            // Check if table entry can be updated with current op (without extend_Fact considered)
            if (hm_table.find(hm_pair)->second <= eval + op_cost) {
            	continue;
            }
        	const int c2 = extend_eval(extend_fact, pre, eval);
        	if (c2 != INT_MAX) {
            	update_hm_entry(hm_pair, c2 + op_cost);
        	}
        }
    }
}

/*
 * Handles changed entries for op. Used if op is achievable with the same cost as in last iteration.
 * In this case, it is enough to consider extending entries with all changed atoms, as potential improvements do not come from op but from entries changed due to oth
 */
void HTwoHeuristic::extend_changed_entry(const OperatorProxy &op) {
    const int op_id = op.get_id();
    const vector<FactPair> &pre = precondition_cache[op_id];
    const int op_base_cost = op.get_cost();
    const int cost = op_cost[op_id];
	for (EffectProxy eff : op.get_effects()) {
    	FactPair effect = eff.get_fact().get_pair();
        for (FactPair entry : changed_entries[op_id]) {
        	Pair hm_pair = effect.var > entry.var ? Pair(entry, effect) : Pair(effect, entry);
            if (hm_table.find(hm_pair)->second <= cost + op_base_cost) {
            	continue;
            }
            int c2 = extend_eval(entry, pre, cost);
        	if (c2 != INT_MAX) {
            	update_hm_entry(hm_pair, c2 + op_base_cost);
        	}
        }
    }
}

/*
 * Evaluates atom set by computing the maximum heuristic value among all its subsets (subset size <= 2). Used for pre(op) and goal.
 */
int HTwoHeuristic::eval(const vector<FactPair> &atom_set) const {
    vector<Pair> pairs = generate_all_pairs(atom_set);
    int max = 0;
    for (Pair &pair : pairs) {
        int h = hm_table.find(pair)->second;
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
 * Evaluates extend_fact + pre. pre already evaluated with integer eval. Runtime is linear in |pre|.
 */
int HTwoHeuristic::extend_eval(const FactPair &extend_fact, const vector<FactPair> &pre, int eval) const {
    int fact_eval = hm_table.find(Pair(extend_fact, FactPair(-1, -1)))->second;
    int max = eval > fact_eval? eval : fact_eval;
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
        int h = hm_table.find(key)->second;

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
 * Adds operators to the queue if atom in pair p is at least partially contained in pre(op). Called if hm_table entry of p is updated.
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
 * Affected operators are potentially added to queue.
 */
void HTwoHeuristic::update_hm_entry(const Pair &p, int val) {
    if (hm_table[p] > val) {
        hm_table[p] = val;
        add_operator_to_queue(p);
    }
}

/*
 * Generates all partial set of size <= 2 from given base atom set.
 */
vector<HTwoHeuristic::Pair> HTwoHeuristic::generate_all_pairs(const vector<FactPair> &base_atom_set) const {
    const size_t n = base_atom_set.size();
	vector<Pair> res;
    res.reserve(n * (n + 1) / 2);

    for (size_t i = 0; i < n; ++i) {
        res.emplace_back(base_atom_set[i], FactPair(-1, -1));

        for (size_t j = i + 1; j < n; ++j) {
            res.emplace_back(base_atom_set[i], base_atom_set[j]);
        }
    }
    return res;
}

void HTwoHeuristic::dump_table() const {
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
