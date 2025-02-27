#include "pi_m_compiled_task.h"

#include "../task_utils/task_properties.h"
#include "../utils/logging.h"


using namespace std;
namespace extra_tasks {

PiMCompiledTask::PiMCompiledTask(const shared_ptr<AbstractTask> &parent) : DelegatingTask(parent) {
    init_meta_atom_map();
    setup_init_and_goal_states();
	setup_meta_operators();

}

void PiMCompiledTask::init_meta_atom_map() {
	meta_atom_map = {{pair(FactPair(-1, -1), FactPair(-1, -1)), 0}};
    int num_var = parent->get_num_variables();
    int index = 1;
    for (int var = 0; var < num_var; ++var) {
        for (int val = 0; val < parent->get_variable_domain_size(var); ++val) {
            for (int value = val; value < parent->get_variable_domain_size(var); ++value) {
                meta_atom_map[pair(FactPair(var, val), FactPair(var, value))] = index;
            	index++;
        	}
            for (int variable = var + 1; variable < num_var; ++variable) {
        		for (int value = 0; value < parent->get_variable_domain_size(variable); ++value) {
                    meta_atom_map[pair(FactPair(var, val), FactPair(variable, value))] = index;
            		index++;
        		}
    		}
        }
    }
}

void PiMCompiledTask::setup_init_and_goal_states() {
    initial_state_values.assign(meta_atom_map.size(), 0);
    vector<int> init_state_values = parent->get_initial_state_values();
    unordered_set<FactPair, FactPairHash> goal_facts;
    for (int i = 0; i < parent->get_num_goals(); ++i) {
    	goal_facts.insert(parent->get_goal_fact(i));
    }

	for (const auto& [current_atoms, i] : meta_atom_map) {
        if (current_atoms.first.var == -1 || init_state_values[current_atoms.first.var] == current_atoms.first.value) {
        	if (current_atoms.second.var == -1 || init_state_values[current_atoms.second.var] == current_atoms.second.value) {
        		initial_state_values[i] = 1;
        	}
        }
        if (current_atoms.first.var == -1 || goal_facts.find(current_atoms.first) != goal_facts.end()) {
        	if (current_atoms.second.var == -1 || goal_facts.find(current_atoms.second) != goal_facts.end()) {
            	goals.push_back(FactPair(i, 1));
        	}
        }
    }
}

void PiMCompiledTask::setup_meta_operators() {
    meta_operators = {};
    for (int op_id = 0; op_id < parent->get_num_operators(); ++op_id) {
        // To check S ∩ (add(o) ∪ del(o)) = ∅
		unordered_set<int> effect_vars;
        vector<FactPair> new_pre = generate_meta_preconditions(op_id);
        vector<FactPair> new_eff = generate_meta_effects(op_id, effect_vars);

        // S = ∅
        meta_operators.push_back(MetaOperator(op_id, FactPair(-1, -1), new_pre, new_eff, parent->get_operator_cost(op_id, false)));

        for (int var = 0; var < parent->get_num_variables(); ++var) {
            if (effect_vars.find(var) != effect_vars.end()) {
            	continue;
            }
            for (int val = 0; val < parent->get_variable_domain_size(var); val++) {
            	FactPair s_atom = FactPair(var, val);
                if (contradict_precondition(op_id, s_atom)) {
                	continue;
                }
            	vector<FactPair> copy_pre(new_pre);
            	copy_pre.push_back(FactPair(meta_atom_map[pair(s_atom, s_atom)], 1));
            	for (int i = 0; i < parent->get_num_operator_preconditions(op_id, false); ++i) {
					FactPair pre = parent->get_operator_precondition(op_id, i, false);
                	FactPair meta_pre = translate_into_meta_atom(pre, s_atom);
                	if (meta_pre.var != -2) {
                    	copy_pre.push_back(meta_pre);
                	}
            	}
            	vector<FactPair> copy_eff(new_eff);
            	for (int i = 0; i < parent->get_num_operator_effects(op_id, false); ++i) {
					FactPair eff = parent->get_operator_effect(op_id, i, false);
                	FactPair meta_eff = translate_into_meta_atom(eff, s_atom);
                	if (meta_eff.var != -2) {
                    	copy_eff.push_back(meta_eff);
                	}
            	}
				meta_operators.push_back(MetaOperator(op_id, s_atom, copy_pre, copy_eff, parent->get_operator_cost(op_id, false)));
            }
        }
    }
}

FactPair PiMCompiledTask::translate_into_meta_atom(FactPair first_atom, FactPair second_atom) {
    pair atom_pair = pair(FactPair(-1, -1), FactPair(-1, -1));
	if (first_atom < second_atom) {
        atom_pair = pair(first_atom, second_atom);
    } else {
    	atom_pair = pair(second_atom, first_atom);
    }
    if (meta_atom_map.find(atom_pair) == meta_atom_map.end()) {
        return FactPair(-2, -2);
    }

    return FactPair(meta_atom_map[atom_pair], 1);
}

bool PiMCompiledTask::contradict_precondition(int op_id, FactPair s_atom) {
	for (int i = 0; i < parent->get_num_operator_preconditions(op_id, false); ++i) {
		FactPair pre = parent->get_operator_precondition(op_id, i, false);
    	if (pre.var == s_atom.var && pre.value != s_atom.value) {
        	return true;
    	}
    }
    return false;
}

vector<FactPair> PiMCompiledTask::generate_meta_preconditions(int op_id) {
    vector<FactPair> new_pre = {FactPair(meta_atom_map[{FactPair(-1, -1), FactPair(-1, -1)}], 1)};
    for (int i = 0; i < parent->get_num_operator_preconditions(op_id, false); ++i) {
		FactPair pre = parent->get_operator_precondition(op_id, i, false);
        for (int j = 0; j < parent->get_num_operator_preconditions(op_id, false); ++j) {
			FactPair second_pre = parent->get_operator_precondition(op_id, j, false);
            FactPair meta_pre = translate_into_meta_atom(pre, second_pre);
            if (meta_pre.var != -2) {
                new_pre.push_back(meta_pre);
            }
        }
    }
    return new_pre;
}

vector<FactPair> PiMCompiledTask::generate_meta_effects(int op_id, unordered_set<int> &effect_vars) {
    vector<FactPair> new_eff;
    for (int i = 0; i < parent->get_num_operator_effects(op_id, false); ++i) {
		FactPair eff = parent->get_operator_effect(op_id, i, false);
        effect_vars.insert(eff.var);
        for (int j = 0; j < parent->get_num_operator_effects(op_id, false); ++j) {
			FactPair second_eff = parent->get_operator_effect(op_id, j, false);
            FactPair meta_eff = translate_into_meta_atom(eff, second_eff);
            if (meta_eff.var != -2) {
                new_eff.push_back(meta_eff);
            }
        }
    }
    return new_eff;
}

int PiMCompiledTask::get_num_variables() const {
    return meta_atom_map.size();
}

int PiMCompiledTask::get_variable_domain_size(int var) const {
	(void)var;
    return 2;
}

int PiMCompiledTask::get_operator_cost(int index, bool is_axiom) const {
    (void)is_axiom;
    return meta_operators[index].cost;
}


int PiMCompiledTask::get_num_operators() const {
    return meta_operators.size();
}

int PiMCompiledTask::get_num_operator_preconditions(int index, bool is_axiom) const {
    (void)is_axiom;
    return meta_operators[index].preconditions.size();
}

FactPair PiMCompiledTask::get_operator_precondition(int op_index, int fact_index, bool is_axiom) const {
    (void)is_axiom;
    return meta_operators[op_index].preconditions[fact_index];
}

int PiMCompiledTask::get_num_operator_effects(int op_index, bool is_axiom) const {
    (void)is_axiom;
	return meta_operators[op_index].effects.size();
}

FactPair PiMCompiledTask::get_operator_effect(int op_index, int eff_index, bool is_axiom) const {
    (void)is_axiom;
    return meta_operators[op_index].effects[eff_index];
}

FactPair PiMCompiledTask::get_goal_fact(int index) const {
    return goals[index];
}

int PiMCompiledTask::get_num_goals() const {
    return goals.size();
}

vector<int> PiMCompiledTask::get_initial_state_values() const {
    return initial_state_values;
}

int PiMCompiledTask::get_num_operator_effect_conditions(
        int op_index, int eff_index, bool is_axiom) const {
    (void) op_index; (void) eff_index; (void)is_axiom;
	return 0;
}



void PiMCompiledTask::convert_state_values_from_parent(std::vector<int> &values) const {
    std::vector<int> new_values(initial_state_values.size(), 0);
	new_values[0] = 1;
	for (size_t i = 0; i < values.size(); ++i) {
		FactPair atom = FactPair(i, values[i]);
		new_values[meta_atom_map.at(pair(atom, atom))] = 1;
		for (size_t j = i + 1; j < values.size(); ++j) {
			FactPair second_atom = FactPair(j, values[j]);
			new_values[meta_atom_map.at(pair(atom, second_atom))] = 1;
		}
	}
	values = std::move(new_values);
}

std::shared_ptr<AbstractTask> build_pi_m_compiled_task(
	const shared_ptr<AbstractTask> &parent) {
	return make_shared<PiMCompiledTask>(parent);
}
}