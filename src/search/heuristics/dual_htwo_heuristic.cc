#include "dual_htwo_heuristic.h"

#include "../plugins/plugin.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"
#include "../tasks/dual_task.h"

#include <cassert>
#include <climits>
#include <set>

using namespace std;

namespace dual_htwo_heuristic {

DualHTwoHeuristic::DualHTwoHeuristic(
    const shared_ptr<AbstractTask> &transform,
    bool cache_estimates, const string &description,
    utils::Verbosity verbosity)
    : HTwoHeuristic(transform, cache_estimates, description, verbosity) {
    dual_task = extra_tasks::build_dual_task(transform);
    log << "Initializied dual task" << endl;
    std::vector<int> values = dual_task->get_initial_state_values();

    TaskProxy original_task_proxy = HTwoHeuristic::task_proxy;
    vector<FactPair> original_goals = HTwoHeuristic::goals;
    // Set task_proxy to dual_task only for creatomg hm_table
    HTwoHeuristic::task_proxy = TaskProxy(*dual_task);
    goals = {};
    // Initialize op caches with dual task (constructor of htwo_heuristic uses original task)
    HTwoHeuristic::init_binary_operators();
    HTwoHeuristic::setup_precondition_of();
    HTwoHeuristic::compute_heuristic(State(*dual_task, std::move(values)));
    HTwoHeuristic::task_proxy = original_task_proxy;
    goals = original_goals;
}


// Removed goal check in the beginning as goals are detected rather fast anyways
int DualHTwoHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = Heuristic::convert_ancestor_state(ancestor_state);
    vector<int> state_values = state.get_unpacked_values();
    dynamic_pointer_cast<extra_tasks::DualTask>(dual_task)->convert_state_values_from_parent(state_values);
	vector<FactPair> state_atoms = {};
    for (size_t i = 0; i < state_values.size(); i++) {
        if (state_values[i] == 1) {
        	state_atoms.push_back(FactPair(i, state_values[i]));
        }
    }
    int h = eval(state_atoms);
    if (h == INT_MAX) {
        return DEAD_END;
    }
    return h;

}

/*
 * Evaluates tuple by computing the maximum heuristic value among all its partial tuples. Used for pre(op) and goal.
 */
int DualHTwoHeuristic::eval(const vector<FactPair> &fp)  {
    vector<Pair> pairs = generate_all_pairs(fp);
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


class DualHTwoHeuristicFeature
    : public plugins::TypedFeature<Evaluator, DualHTwoHeuristic> {
public:
    DualHTwoHeuristicFeature() : TypedFeature("dualh2") {
        document_title("Dual h^2 heuristic");

        add_heuristic_options_to_feature(*this, "dualh2");

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

    virtual shared_ptr<DualHTwoHeuristic> create_component(
        const plugins::Options &opts,
        const utils::Context &) const override {
        return plugins::make_shared_from_arg_tuples<DualHTwoHeuristic>(
            get_heuristic_arguments_from_options(opts)
            );
    }
};

static plugins::FeaturePlugin<DualHTwoHeuristicFeature> _plugin;
}
