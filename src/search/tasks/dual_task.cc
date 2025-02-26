#include "dual_task.h"

#include "../task_utils/task_properties.h"

using namespace std;
namespace extra_tasks {

DualTask::DualTask(const std::shared_ptr<AbstractTask> &parent) : DelegatingTask(parent) {
	init_strips_variables();
    setup_init_and_goal_states();
    setup_dual_operators();

    print_task();
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
	dual_operator_pre = {};
    dual_operator_eff = {};
    for (int op_id = 0; op_id < parent->get_num_operators(); ++op_id) {
        dual_operator_pre.push_back({});
        dual_operator_eff.push_back({});
        for (int var = 0; var < parent->get_num_operator_preconditions(op_id, false); ++var) {
			FactPair pre = parent->get_operator_precondition(op_id, var, false);
			for (int val = 0; val < parent->get_variable_domain_size(pre.var); ++val) {
                // Add precondition contradictions to effects
				if (pre.value != val) {
                	dual_operator_eff[op_id].push_back(FactPair(strips_atom_map[FactPair(pre.var, val)], 1));
                }
			}
        }
    	for (int var = 0; var < parent->get_num_operator_effects(op_id, false); ++var) {
			FactPair eff = parent->get_operator_effect(op_id, var, false);
            FactPair new_eff = FactPair(strips_atom_map[eff], 1);
            if (find(dual_operator_eff[op_id].begin(), dual_operator_eff[op_id].end(), new_eff) == dual_operator_eff[op_id].end()) {
            	dual_operator_eff[op_id].push_back(new_eff);
            }
			for (int val = 0; val < parent->get_variable_domain_size(eff.var); ++val) {
                // Add delete effects (atoms that contradict eff) to preconditions
				if (eff.value != val) {
                	dual_operator_pre[op_id].push_back(FactPair(strips_atom_map[FactPair(eff.var, val)], 1));
                }
			}
        }
    }
}

void DualTask::print_task() {
	cout << "Task " << endl;

	for (int i = 0; i < parent->get_num_operators(); ++i) {
    	cout << "Operator " << i << endl << "pre: ";
        for (int var = 0; var < parent->get_num_operator_preconditions(i, false); ++var) {
        	cout << parent->get_operator_precondition(i, var, false) << ", ";
        }
        cout << endl << "eff: ";
        for (int var = 0; var < parent->get_num_operator_effects(i, false); ++var) {
        	cout << parent->get_operator_effect(i, var, false) << ", ";
        }
        cout << endl;
    }
    cout << endl << "Initial state: " << endl;
    for (int i = 0; i < parent->get_num_variables(); ++i) {
    	cout << i << "=" << parent->get_initial_state_values()[i] << ", ";
    }
    cout << endl << "Goal: " << endl;
    for (int i = 0; i < parent->get_num_goals(); ++i) {
    	cout << parent->get_goal_fact(i) << ", ";
    }

    cout << endl << endl << "Dual Task: " << endl;
	for (const auto& [atom, new_var] : strips_atom_map) {
    	cout << atom << " -> " << new_var << endl;
    }
	cout << endl;
    for (int i = 0; i < dual_operator_pre.size(); ++i) {
    	cout << "Operator " << i << endl << "pre: ";
        for (int var = 0; var < dual_operator_pre[i].size(); ++var) {
        	cout << dual_operator_pre[i][var] << ", ";
        }
        cout << endl << "eff: ";
        for (int var = 0; var < dual_operator_eff[i].size(); ++var) {
        	cout << dual_operator_eff[i][var] << ", ";
        }
        cout << endl;
    }
    cout << endl << "Initial state: " << endl;
    for (int i = 0; i < initial_state_values.size(); ++i) {
    	cout << i << "=" << initial_state_values[i] << ", ";
    }
    cout << endl << "Goals: " << endl;
    for (int i = 0; i < goals.size(); ++i) {
    	cout << goals[i] << ", ";
    }
    cout << endl;
}


int DualTask::get_num_variables() const {
    return strips_atom_map.size();
}

int DualTask::get_variable_domain_size(int var) const {
	(void)var;
    return 2;
}


int DualTask::get_num_operators() const {
    return dual_operator_pre.size();
}

int DualTask::get_num_operator_preconditions(int index, bool is_axiom) const {
    (void)is_axiom;
    return dual_operator_pre[index].size();
}

FactPair DualTask::get_operator_precondition(int op_index, int fact_index, bool is_axiom) const {
    (void)is_axiom;
    return dual_operator_pre[op_index][fact_index];
}

int DualTask::get_num_operator_effects(int op_index, bool is_axiom) const {
    (void)is_axiom;
	return dual_operator_eff[op_index].size();
}

FactPair DualTask::get_operator_effect(int op_index, int eff_index, bool is_axiom) const {
    (void)is_axiom;
    return dual_operator_eff[op_index][eff_index];
}

FactPair DualTask::get_goal_fact(int index) const {
    return goals[index];
}

int DualTask::get_num_goals() const {
    return goals.size();
}

vector<int> DualTask::get_initial_state_values() const {
    return initial_state_values;
}

int DualTask::get_num_operator_effect_conditions(
        int op_index, int eff_index, bool is_axiom) const {
    (void) op_index; (void) eff_index; (void)is_axiom;
	return 0;
}

void DualTask::convert_state_values_from_parent(std::vector<int> &values) const {
	for (size_t i = 0; i < strips_atom_map.size(); ++i) {
		if (values[i] == 0) {
        	values[i] = 1;
        } else {
            values[i] = 0;
        }
	}
}

std::shared_ptr<AbstractTask> build_dual_task(
	const shared_ptr<AbstractTask> &parent) {
	return make_shared<DualTask>(parent);
}
}