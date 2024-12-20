#include "hm_heuristic.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"

#include <cassert>
#include <limits>
#include <set>

using namespace std;

namespace hm_heuristic {

/*
 * Constructor for the HMHeuristic class.
 * Precomputes all possible tuples of size <= m.
 */
HMHeuristic::HMHeuristic(
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
    generate_all_tuples();
}


bool HMHeuristic::dead_ends_are_reliable() const {
    return !task_properties::has_axioms(task_proxy) && !has_cond_effects;
}


/*
 * Computes the h^m value for a given state:
 * Checks if state is a goal state (heuristic = 0 if true). Initializes h^m table with state facts.
 * Updates h^m table to propagate values. Evaluates goal facts to compute the heuristic value.
 */
int HMHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    } else {
        Tuple s_tup = task_properties::get_fact_pairs(state);

        init_hm_table(s_tup);
        update_hm_table();

        int h = eval(goals);

        if (h == numeric_limits<int>::max())
            return DEAD_END;
        return h;
    }
}

/*
 * Initializes h^m table.
 * If tuple is contained in input tuple assigns 0, and infinity otherwise.
 */
void HMHeuristic::init_hm_table(const Tuple &t) {
    for (auto &hm_ent : hm_table) {
        const Tuple &tuple = hm_ent.first;
        int h_val = check_tuple_in_tuple(tuple, t);
        hm_table[tuple] = h_val;
    }
}

/*
 * Iteratively updates the h^m table until no further improvements are made.
 * For each operator:
 * - Computes the h value of preconditions.
 * - If preconditions are reachable generates all partial effect tuples for table update.
 * - Extends partial tuples with size < m to explore potential improvements.
 * - Runtime: O(|O|mV^4m log V) with following assumptions:
 *  - |pre(o)|, |eff(o)| = O(V)
 *  - O(V^m) updates of hm_table
 */
void HMHeuristic::update_hm_table() { // O(|O|mV^(4m) log V)
    do {  // O(|O|mV^4m log V) if num iterations = O(V^m)
        //dump_table();
        was_updated = false;
        //log << endl << "New Loop" << endl;
        for (OperatorProxy op : task_proxy.get_operators()) { // O(|O|mV^3m log V)
            Tuple pre = get_operator_pre(op);
            //log << "Operator " << op.get_id() << " pre: " << get_operator_pre(op) << " eff: " << get_operator_eff(op) << endl;

            int c1 = eval(pre);
            if (c1 != numeric_limits<int>::max()) {
                Tuple eff = get_operator_eff(op);
                vector<Tuple> partial_effs;
                generate_all_partial_tuples(eff, partial_effs); // O(V^m) if |eff| = O(V)
                for (Tuple &partial_eff : partial_effs) { // O(mV^3m log V)
                    //if (c1 + op.get_cost() < hm_table[partial_eff]) {
                        //if (partial_eff.size() == 1) {
                           //log << "Eff Update: ([" << partial_eff[0].var << "=" << partial_eff[0].value << "]) = " << c1 << " + " << op.get_cost() << endl;
                        //} else {
                            //log << "Eff Update: ([" << partial_eff[0].var << "=" << partial_eff[0].value << "," << partial_eff[1].var << "=" << partial_eff[1].value << "]) = " << c1 << " + " << op.get_cost() << endl;
                        //}
                    //}
                    update_hm_entry(partial_eff, c1 + op.get_cost());

                    int eff_size = partial_eff.size();
                    if (eff_size < m) {
                        extend_tuple(partial_eff, op); // O(mV^2m log V)
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
void HMHeuristic::extend_tuple(const Tuple &t, const OperatorProxy &op) { // O(mV^2m log V)
     //log << "Extend Tuple " << t[0].var << "=" << t[0].value << endl;
    for (auto &hm_ent : hm_table) { // O(mV^2m log V)
        const Tuple &tuple = hm_ent.first;
        bool contradict = false;
        for (const FactPair &fact : tuple) { // O(mV log V)
            if (contradict_effect_of(op, fact.var, fact.value)) { // O(V log V)
                contradict = true;
                break;
            }
        }
        // only if t is fully contained in tuple
        if (!contradict && (tuple.size() > t.size()) && (check_tuple_in_tuple(t, tuple) == 0)) { // V = {a, b, c}, o = {a, b}, hm_table(b, c) -> pre = {a, c)
            Tuple pre = get_operator_pre(op); // O(V log V)

            Tuple others;
            for (const FactPair &fact : tuple) { // O(V log V)
                // if fact not contained in t
                if (find(t.begin(), t.end(), fact) == t.end()) {
                    others.push_back(fact);
                    // if fact not contained in pre
                    if (find(pre.begin(), pre.end(), fact) == pre.end()) {
                        pre.push_back(fact);
                    }
                }
            }

            sort(pre.begin(), pre.end()); // O(V log V)

            // Checks if no duplicate fact var in pre
            set<int> vars;
            bool is_valid = true;
            for (const FactPair &fact : pre) { // O(V log V)
                if (vars.count(fact.var) != 0) { // O(log V)
                    is_valid = false;
                    break;
                }
                vars.insert(fact.var);
            }
            // Update Table
            if (is_valid) {
                int c2 = eval(pre); // O(mV^m log V) if |pre| = O(V)
                if (c2 != numeric_limits<int>::max()) {
                    //if (c2 + op.get_cost() < hm_table[tuple]) {
                        //if (tuple.size() == 1) {
                            //log << "Ext Update: ([" << tuple[0].var << "=" << tuple[0].value << "]) = " << c2 << " + " << op.get_cost() << endl;
                        //} else {
                            //log << "Ext Update: ([" << tuple[0].var << "=" << tuple[0].value << "," << tuple[1].var << "=" << tuple[1].value << "]) = " << c2 << " + " << op.get_cost() << endl;
                        //}
                    //}
                    update_hm_entry(tuple, c2 + op.get_cost()); // O(m log V)
                }
            }
        }
    }
}


/*
 * Evaluates tuple by computing the maximum heuristic value among all its partial tuples.
 */
int HMHeuristic::eval(const Tuple &t) const {  // O(mV^m log V) 
    vector<Tuple> partial;
    generate_all_partial_tuples(t, partial); // O(V^m)
    int max = 0;
    for (Tuple &tuple : partial) { // O (mV^m log V)
        assert(hm_table.count(tuple) == 1); // O(log V^m)

        int h = hm_table.at(tuple); // O(log V^m)
        if (h > max) {
            max = h;
        }
    }
    return max;
}

/*
 * Updates the heuristic value of a tuple in the h^m table.
 * Sets "was_updated" flag to true to indicate a change occurred.
 */
int HMHeuristic::update_hm_entry(const Tuple &t, int val) {
    assert(hm_table.count(t) == 1);
    if (hm_table[t] > val) {
        hm_table[t] = val;
        was_updated = true;
    }
    return val;
}


/*
 * Checks if tuple is fully contained in another tuple.
 * Returns 0 if fully contained, and infinity otherwise.
 */
int HMHeuristic::check_tuple_in_tuple(
    const Tuple &tuple, const Tuple &big_tuple) const {
    for (const FactPair &fact0 : tuple) {
        bool found = false;
        for (auto &fact1 : big_tuple) {
            if (fact0 == fact1) {
                found = true;
                break;
            }
        }
        if (!found) {
            return numeric_limits<int>::max();
        }
    }
    return 0;
}

HMHeuristic::Tuple HMHeuristic::get_operator_pre(const OperatorProxy &op) const {
    Tuple preconditions = task_properties::get_fact_pairs(op.get_preconditions());
    sort(preconditions.begin(), preconditions.end());
    return preconditions;
}


HMHeuristic::Tuple HMHeuristic::get_operator_eff(const OperatorProxy &op) const {
    Tuple effects;
    for (EffectProxy eff : op.get_effects()) {
        effects.push_back(eff.get_fact().get_pair());
    }
    sort(effects.begin(), effects.end());
    return effects;
}


bool HMHeuristic::contradict_effect_of(
    const OperatorProxy &op, int var, int val) const {
    for (EffectProxy eff : op.get_effects()) {
        FactProxy fact = eff.get_fact();
        if (fact.get_variable().get_id() == var && fact.get_value() != val) {
            return true;
        }
    }
    return false;
}

/*
 * Recursively generates all possible tuples of size <= m over variables of the task.
 * All possible states with <= m variables are stored in hm_table.
 */
void HMHeuristic::generate_all_tuples() {
    Tuple t;
    generate_all_tuples_aux(0, m, t);
}


void HMHeuristic::generate_all_tuples_aux(int var, int sz, const Tuple &base) {
    int num_variables = task_proxy.get_variables().size();
    for (int i = var; i < num_variables; ++i) {
        int domain_size = task_proxy.get_variables()[i].get_domain_size();
        for (int j = 0; j < domain_size; ++j) {
            Tuple tuple(base);
            tuple.emplace_back(i, j);
            hm_table[tuple] = 0;
            if (sz > 1) {
                generate_all_tuples_aux(i + 1, sz - 1, tuple);
            }
        }
    }
}


/*
 * Generates all partial tuples of size <= m from given base tuple.
 */
void HMHeuristic::generate_all_partial_tuples(
    const Tuple &base_tuple, vector<Tuple> &res) const {
    Tuple t;
    generate_all_partial_tuples_aux(base_tuple, t, 0, m, res);
}


void HMHeuristic::generate_all_partial_tuples_aux(
    const Tuple &base_tuple, const Tuple &t, int index, int sz, vector<Tuple> &res) const {
    if (sz == 1) {
        for (size_t i = index; i < base_tuple.size(); ++i) {
            Tuple tuple(t);
            tuple.push_back(base_tuple[i]);
            res.push_back(tuple);
        }
    } else {
        for (size_t i = index; i < base_tuple.size(); ++i) {
            Tuple tuple(t);
            tuple.push_back(base_tuple[i]);
            res.push_back(tuple);
            generate_all_partial_tuples_aux(base_tuple, tuple, i + 1, sz - 1, res);
        }
    }
}


void HMHeuristic::dump_table() const {
  	stringstream ss;
    for (auto &hm_ent : hm_table) {
        if (hm_ent.second == numeric_limits<int>::max()) {
          	continue;
      	}
        ss << "h(" << hm_ent.first << ") = " << hm_ent.second << endl;
    }
    log << ss.str() << endl;
}

class HMHeuristicFeature
    : public plugins::TypedFeature<Evaluator, HMHeuristic> {
public:
    HMHeuristicFeature() : TypedFeature("hm") {
        document_title("h^m heuristic");

        add_option<int>("m", "subset size", "2", plugins::Bounds("1", "infinity"));
        add_heuristic_options_to_feature(*this, "hm");

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

    virtual shared_ptr<HMHeuristic> create_component(
        const plugins::Options &opts,
        const utils::Context &) const override {
        return plugins::make_shared_from_arg_tuples<HMHeuristic>(
            opts.get<int>("m"),
            get_heuristic_arguments_from_options(opts)
            );
    }
};

static plugins::FeaturePlugin<HMHeuristicFeature> _plugin;
}
