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
    TaskProxy original_task_proxy = task_proxy;
    task_proxy = TaskProxy(*dual_task);
    std::vector<int> values = dual_task->get_initial_state_values();
	log << endl;
    HTwoHeuristic::compute_heuristic(State(*dual_task, std::move(values)));
    dual_hm_table = HTwoHeuristic::hm_table;
	//HTwoHeuristic::print_table();
    task_proxy = original_task_proxy;
}



int DualHTwoHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    if (task_properties::is_goal_state(task_proxy, state)) {
        return 0;
    }
    return 0;

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
