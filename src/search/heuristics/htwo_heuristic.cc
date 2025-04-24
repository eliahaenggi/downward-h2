#include "htwo_heuristic.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cassert>
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
    init_binary_operators();
    setup_precondition_of();
    vector<Pair> par_goals = generate_all_pairs(goals);
    // Used for checks if atom pair is part of goal
    partial_goals = unordered_set<Pair, PairHash>(par_goals.begin(), par_goals.end());
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
    int h = update_hm_table();
    if (h == INT_MAX) {
        return DEAD_END;
    }
    return h;
}

/**
* Sets up all auxiliary data structures concerning operators.
*/
void HTwoHeuristic::init_binary_operators() {
	binary_operators = {};

	for (OperatorProxy op : task_proxy.get_operators()) {

        vector<FactPair> preconditions = task_properties::get_fact_pairs(op.get_preconditions());
    	sort(preconditions.begin(), preconditions.end());

    	vector<FactPair> effects;
        // Used for fast checks if var is an effect of op
        vector<bool> effect_conflict = std::vector<bool>(task_proxy.get_variables().size(), false);
    	for (EffectProxy eff : op.get_effects()) {
        	effects.push_back(eff.get_fact().get_pair());
            effect_conflict[eff.get_fact().get_pair().var] = true;
    	}
    	sort(effects.begin(), effects.end());
        std::vector<Pair> partial_effects = generate_all_pairs(effects);
		//log << "Op " << op.get_id() << ", pre: " << preconditions << ", eff: " << effects << endl;
        for (auto& partial : partial_effects) {
            //log << "Add bin op pre: " << preconditions << ", eff: " << partial.first << ", " << partial.second << endl;
        	binary_operators.push_back(BinaryOperator(preconditions, partial, op.get_cost(), op.get_id(), binary_operators.size()));

            if (partial.second.var == -1) {
            	// Extend binary operators with preconditions still true
                for (auto pre : preconditions) {
                	if (!effect_conflict[pre.var]) {
						Pair pair = partial.first.var < pre.var ? Pair(partial.first, pre) : Pair(pre, partial.first);
                        //log << "Add bin op pre: " << preconditions << ", eff: " << pair.first << ", " << pair.second << endl;
        				binary_operators.push_back(BinaryOperator(preconditions, pair, op.get_cost(), op.get_id(), binary_operators.size()));
                	}
                }
                // Extend binary operators with other atoms
                extend_binary_operators(partial.first, op, effect_conflict);
            }
        }
    }
}

void HTwoHeuristic::extend_binary_operators(const FactPair &f, const OperatorProxy &op, vector<bool>& effect_conflict) {
    const auto &variables = task_proxy.get_variables();
    const int op_id = op.get_id();
    const int op_cost = op.get_cost();
    for (size_t i = 0; i < variables.size(); ++i) {
        if (effect_conflict[i]) {
        	continue;
        }
        const int domain_size = variables[i].get_domain_size();
    	for (int j = 0; j < domain_size; ++j) {
        	const FactPair extend_fact = FactPair(i, j);
            const Pair hm_pair = f.var > extend_fact.var ? Pair(extend_fact, f) : Pair(f, extend_fact);
            vector<FactPair> pre = task_properties::get_fact_pairs(op.get_preconditions());
            if (find(pre.begin(), pre.end(), extend_fact) == pre.end()) {
            	pre.push_back(extend_fact);
            }
            unordered_set<int> vars;
            bool is_valid = true;
            for (const FactPair &fact : pre) {
                if (vars.count(fact.var) != 0) {
                    is_valid = false;
                    break;
                }
                vars.insert(fact.var);
            }
            sort(pre.begin(), pre.end());
            if (is_valid) {
				binary_operators.push_back(BinaryOperator(pre, hm_pair, op_cost, op_id, binary_operators.size()));
            }
        }
    }
}

void HTwoHeuristic::setup_precondition_of() {
    const int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        const int domain1_size = task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < domain1_size; ++j) {
            Pair single_pair(FactPair(i, j), FactPair(-1, -1));
            precondition_of[single_pair] = {};
            for (auto& op : binary_operators) {
            	if (find(op.preconditions.begin(), op.preconditions.end(), FactPair(i, j)) != op.preconditions.end()) {
                	precondition_of[single_pair].push_back(op.id);
                }
            }
            for (int k = i + 1; k < num_variables; ++k) {
                const int domain2_size = task_proxy.get_variables()[k].get_domain_size();
                for (int l = 0; l < domain2_size; ++l) {
                    Pair pair(FactPair(i, j), FactPair(k, l));
                    precondition_of[pair] = {};
                    for (auto& op : binary_operators) {
                    	if (find(op.preconditions.begin(), op.preconditions.end(), pair.first) != op.preconditions.end()) {
                            if (find(op.preconditions.begin(), op.preconditions.end(), pair.second) != op.preconditions.end()) {
                        		precondition_of[pair].push_back(op.id);
                    		}
                    	}
                    }
                }
            }
        }
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
    queue.clear();

    const int num_variables = task_proxy.get_variables().size();
    for (int i = 0; i < num_variables; ++i) {
        const int domain1_size = task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < domain1_size; ++j) {
            Pair single_pair(FactPair(i, j), FactPair(-1, -1));
            int init_h = check_in_initial_state(single_pair, state_facts_set);
            hm_table[single_pair] = INT_MAX;
            // If part of init state
            if (init_h == 0) {
            	enqueue_if_necessary(single_pair, init_h);
            }

            for (int k = i + 1; k < num_variables; ++k) {
                const int domain2_size = task_proxy.get_variables()[k].get_domain_size();
                for (int l = 0; l < domain2_size; ++l) {
                    Pair pair(FactPair(i, j), FactPair(k, l));
            		init_h = check_in_initial_state(pair, state_facts_set);
            		hm_table[pair] = INT_MAX;
                    // If part of init state
                    if (init_h == 0) {
                    	enqueue_if_necessary(pair, init_h);
                    }
                }
            }
        }
    }
	// Set unsatisfied pre and cost of binary operators, insert precondition free operators to queue
    for (BinaryOperator &op : binary_operators) {
    	op.unsatisfied_preconditions = generate_all_pairs(op.preconditions).size();
        //log << "Op " << op.old_op_id << " : " << op.id << " " << op.preconditions << endl;
        op.cost = op.base_cost;
        if (op.unsatisfied_preconditions == 0) {
        	enqueue_if_necessary(op.effect, op.base_cost);
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


int HTwoHeuristic::update_hm_table() {
    int unsolved_goals = generate_all_pairs(goals).size();
    while (!queue.empty()) {
      	pair<int, Pair> top_pair = queue.pop();
        int distance = top_pair.first;
        Pair pair = top_pair.second;
        //log << pair.first << ", " << pair.second << "   " << distance << endl;
        if (hm_table.at(pair) < distance) {
        	continue;
        }
        if (partial_goals.find(pair) != partial_goals.end() && --unsolved_goals == 0) {
            return distance;
        }
        for (auto& op_id : precondition_of[pair]) {
            BinaryOperator &op = binary_operators[op_id];
            --op.unsatisfied_preconditions;
            if (op.unsatisfied_preconditions == 0) {
                op.cost = max(op.cost, hm_table.at(pair) + op.base_cost);
                enqueue_if_necessary(op.effect, op.cost);
            }
        }
    }
    return INT_MAX;
}

/*
 * Generates all partial set of size <= 2 from given base tuple.
 */
vector<HTwoHeuristic::Pair> HTwoHeuristic::generate_all_pairs(const vector<FactPair> &base_tuple) const {
    const size_t n = base_tuple.size();
	vector<Pair> res;
    res.reserve(n * (n + 1) / 2);

    for (size_t i = 0; i < n; ++i) {
        res.emplace_back(base_tuple[i], FactPair(-1, -1));

        for (size_t j = i + 1; j < n; ++j) {
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
