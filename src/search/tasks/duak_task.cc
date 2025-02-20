#include "dual_task.h"

#include "../task_utils/task_properties.h"

using namespace std;
namespace extra_tasks {

DualTask::DualTask(const std::shared_ptr<AbstractTask> &parent) : DelegatingTask(parent) {
	init_strips_variables();
    setup_init_and_goal_states();
    setup_dual_operators();
}

void DualTask::init_strips_variables() {
	strips_atom_map.clear();
    int index = 0;
    for (int var = 0; var < parent->get_num_variables(); ++var) {
        for (int val = 0; val < parent->get_variable_domain_size(var); ++val) {
        	strips_atom_map[FactPair(var, val)] = index;
            index++;
        }
    }
}

void DualTask::setup_init_and_goal_states() {
    initial_state_values.assign(strips_atom_map.size(), 0);
    vector<int> init_state_values = parent->get_initial_state_values();
    unordered_set<FactPair, FactPairHash> goal_facts;
    for (int i = 0; i < parent->get_num_goals(); ++i) {
    	goal_facts.insert(parent->get_goal_fact(i));
    }

	for (const auto& [atom, i] : strips_atom_map) {
        if (init_state_values[atom.var] != atom.value) {
        	goals.push_back(FactPair(i, 1));
        }
        if (goal_facts.find(atom) == goal_facts.end()) {
        	initial_state_values[i] = 1;
        }
    }
}

void DualTask::setup_dual_operators() {
	dual_operator_pre.clear();
    dual_operator_eff.clear();
    for (int op_id = 0; op_id < parent->get_num_operators(); ++op_id) {
        dual_operator_pre[op_id] = {};
        dual_operator_eff[op_id] = {};
    	for (int var = 0; var < parent->get_num_operator_effects(op_id, false); ++var) {
			FactPair eff = parent->get_operator_effect(op_id, var, false);
            dual_operator_eff[op_id].push_back(eff);
			for (int val = 0; val < parent->get_variable_domain_size(eff.var); ++val) {
                // Add delete effects (atoms that contradict eff) to preconditions
				if (eff.value != val) {
                	dual_operator_pre[op_id].push_back(eff);
                }
			}
        }

        for (int var = 0; var < parent->get_num_operator_preconditions(op_id, false); ++var) {
			FactPair pre = parent->get_operator_precondition(op_id, var, false);
            dual_operator_eff[op_id].push_back(pre);
			for (int val = 0; val < parent->get_variable_domain_size(pre.var); ++val) {
                // Add precondition contradictions to effects
				if (pre.value != val) {
                	dual_operator_eff[op_id].push_back(pre);
                }
			}
        }
    }
}

std::shared_ptr<AbstractTask> build_dual_task(
	const shared_ptr<AbstractTask> &parent) {
	return make_shared<DualTask>(parent);
}
}