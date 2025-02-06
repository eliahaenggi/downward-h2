#include "pi_m_compiled_task.h"

#include "../task_utils/task_properties.h"
#include <unordered_set>



using namespace std;
namespace extra_tasks {

PiMCompiledTask::PiMCompiledTask(const shared_ptr<AbstractTask> &parent) : DelegatingTask(parent) {
  	TaskProxy parent_proxy(*parent);
    if (parent->get_num_axioms() > 0) {
        ABORT("DomainAbstractedTask doesn't support axioms.");
    }
	if (task_properties::has_conditional_effects(parent_proxy)) {
		ABORT("Pi^m Compilation doesn't support conditional effects.");
	}

    store_old_ops();
    init_meta_atom_map();
    setup_init_and_goal_states();
	setup_new_ops();

    dump_compiled_task();
}

void PiMCompiledTask::store_old_ops() {
	int num_op = parent->get_num_operators();
    for (int op_id = 0; op_id < num_op; ++op_id) {
    	vector<FactPair> old_pre_op;
        for (int pre_id = 0; pre_id < parent->get_num_operator_preconditions(op_id, false); pre_id++) {
        	old_pre_op.push_back(parent->get_operator_precondition(op_id, pre_id, false));
        }
        old_pre.push_back(old_pre_op);
      	vector<FactPair> old_eff_op;
        for (int eff_id = 0; eff_id < parent->get_num_operator_effects(op_id, false); eff_id++) {
        	old_eff_op.push_back(parent->get_operator_effect(op_id, eff_id, false));
        }
        old_eff.push_back(old_eff_op);
    }
}

void PiMCompiledTask::init_meta_atom_map() {
	meta_atom_map = {{pair(FactPair(-1, -1), FactPair(-1, -1)), 0}};
    fact_names = {{"not v_∅", "v_∅"}};
    int num_var = parent->get_num_variables();
    int index = 1;
    for (int var = 0; var < num_var; ++var) {
        for (int val = 0; val < parent->get_variable_domain_size(var); ++val) {
            string var_name = "v_" + to_string(var) + "=" + to_string(val);
            for (int variable = var; variable < num_var; ++variable) {
        		for (int value = 0; value < parent->get_variable_domain_size(variable); ++value) {
                    meta_atom_map[pair(FactPair(var, val), FactPair(variable, value))] = index;
            		index++;
                    if (FactPair(var, val) ==  FactPair(variable, value)) {
                    	fact_names.push_back({"not " + var_name, var_name});
                    } else {
                    	string full_var_name = var_name + "," + to_string(variable) + "=" + to_string(value);
                    	fact_names.push_back({"not " + full_var_name, full_var_name});
                    }
        		}
    		}
        }
    }
}

void PiMCompiledTask::setup_init_and_goal_states() {
	domain_size.assign(meta_atom_map.size(), 2);
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

void PiMCompiledTask::setup_new_ops() {
    op_pre = {};
    op_eff = {};
    op_cost = {};
    op_list = {};
    for (int op_id = 0; op_id < parent->get_num_operators(); ++op_id) {
        // Handle S = ∅
        op_list.push_back(pair(op_id, FactPair(-1, -1)));
        op_cost.push_back(parent->get_operator_cost(op_id, false));
        vector<FactPair> old_op_preconditions = old_pre[op_id];
        vector<FactPair> new_pre = {FactPair(meta_atom_map[pair(FactPair(-1,-1), FactPair(-1,-1))], 1)};
		for (size_t pre_id = 0; pre_id < old_op_preconditions.size(); ++pre_id) {
        	FactPair pre = old_op_preconditions[pre_id];
            for (size_t second_pre_id = pre_id; second_pre_id < old_op_preconditions.size(); ++second_pre_id) {
                FactPair second_pre = old_op_preconditions[second_pre_id];
                FactPair meta_pre = translate_into_meta_atom(pre, second_pre);
                if (meta_pre.var != -2) {
                    new_pre.push_back(meta_pre);
                }
            }
		}
        op_pre.push_back(new_pre);
        // To Check S ∩ (add(o) ∪ del(o)) = ∅
		unordered_set<int> effect_vars;
       	vector<FactPair> old_op_effects = old_eff[op_id];
        vector<FactPair> new_eff = {};
        for (size_t eff_id = 0; eff_id < old_op_effects.size(); ++eff_id) {
            FactPair eff = old_op_effects[eff_id];
        	effect_vars.insert(eff.var);
            for (size_t second_eff_id = eff_id; second_eff_id < old_op_effects.size(); ++second_eff_id) {
            	FactPair second_eff = old_op_effects[second_eff_id];
                FactPair meta_eff = translate_into_meta_atom(eff, second_eff);
                if (meta_eff.var != -2) {
                    new_eff.push_back(meta_eff);
                }
            }
        }
        op_eff.push_back(new_eff);
        for (int var = 0; var < parent->get_num_variables(); ++var) {
            // S ∩ (add(o) ∪ del(o)) = ∅
            if (effect_vars.find(var) != effect_vars.end()) {
            	continue;
            }
            for (int val = 0; val < 2; val++) {
            	FactPair s_atom = FactPair(var, val);
                /*
            	if (contradict_precondition(op_id, s_atom)) {
            		continue;
            	}
                 */
            	op_list.push_back(pair(op_id, s_atom));
            	op_cost.push_back(parent->get_operator_cost(op_id, false));
            	vector<FactPair> copy_pre(new_pre);
            	copy_pre.push_back(FactPair(meta_atom_map[pair(s_atom, s_atom)], 1));
            	for (FactPair pre : old_op_preconditions) {
                	FactPair meta_pre = translate_into_meta_atom(pre, s_atom);
                	if (meta_pre.var != -2) {
                    	copy_pre.push_back(meta_pre);
                	}
            	}
            	op_pre.push_back(copy_pre);
            	vector<FactPair> copy_eff(new_eff);
        		for (FactPair eff : old_op_effects) {
                	FactPair meta_eff = translate_into_meta_atom(eff, s_atom);
                	if (meta_eff.var != -2) {
                    	copy_eff.push_back(meta_eff);
                	}
            	}
				op_eff.push_back(copy_eff);
            }
        }
    }
}

FactPair PiMCompiledTask::translate_into_meta_atom(FactPair first_atom, FactPair second_atom) {
    pair atom_pair = pair(FactPair(-1, -1), FactPair(-1, -1));
	if (first_atom.var < second_atom.var) {
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
	for (FactPair pre : old_pre[op_id]) {
    	if (pre.var == s_atom.var && pre.value != s_atom.value) {
        	return true;
    	}
    }
    return false;
}

void PiMCompiledTask::dump_compiled_task() {
      for (size_t i = 0; i < op_list.size(); ++i) {
    	cout << get_operator_name(i, false) << ", pre: ";
        for (auto pre : old_pre[op_list[i].first]) {
        	cout << pre << ", ";
        }
        cout << " eff: ";
        for (auto eff : old_eff[op_list[i].first]) {
        	cout << eff << ", ";
        }
        cout << endl << "pre: ";
        for (auto pre : op_pre[i]) {
        	cout << fact_names[pre.var][pre.value] <<  " , ";
        }
        cout << endl << "eff: ";
        for (auto eff : op_eff[i]) {
        	cout << fact_names[eff.var][eff.value] << " , ";
        }
        cout << endl << endl;
    }
    cout << endl << "Init state: ";
    for (size_t i = 0; i < parent->get_initial_state_values().size(); i++) {
    	cout << i << "=" << parent->get_initial_state_values()[i] << ", ";
    }
    cout << endl << "Compiled init state: ";
    for (size_t i = 0; i < initial_state_values.size(); i++) {
        if (initial_state_values[i] == 1) {
    		cout << fact_names[i][1] << ", ";
        }
    }
    cout << endl << endl << "Goals: ";
    for (int i = 0; i < parent->get_num_goals(); i++) {
    	cout << parent->get_goal_fact(i) << ", ";
    }
    cout << endl << "Compiled goal state: ";
    for (size_t i = 0; i < goals.size(); i++) {
    		cout << fact_names[goals[i].var][1] << ", ";
    }
    cout << endl << endl;
}


int PiMCompiledTask::get_num_variables() const {
    return domain_size.size();
}

string PiMCompiledTask::get_variable_name(int var) const {
    return fact_names[var][1];
}

int PiMCompiledTask::get_variable_domain_size(int var) const {
    return domain_size[var];
}

int PiMCompiledTask::get_operator_cost(int index, bool is_axiom) const {
    (void)is_axiom;
    return op_cost[index];
}

string PiMCompiledTask::get_fact_name(const FactPair &fact) const {
    return fact_names[fact.var][fact.value];
}

string PiMCompiledTask::get_operator_name(int index, bool is_axiom) const {
    (void)is_axiom;
    if (op_list[index].second.var == -1) {
    	return "o_" + to_string(op_list[index].first) + ",∅";
    }
    return "o_" + to_string(op_list[index].first) + "," +
           to_string(op_list[index].second.var) + "=" +
           to_string(op_list[index].second.value);
}

int PiMCompiledTask::get_num_operators() const {
    return op_list.size();
}

int PiMCompiledTask::get_num_operator_preconditions(int index, bool is_axiom) const {
    (void)is_axiom;
    return op_pre[index].size();
}

FactPair PiMCompiledTask::get_operator_precondition(int op_index, int fact_index, bool is_axiom) const {
    (void)is_axiom;
    return op_pre[op_index][fact_index];
}

int PiMCompiledTask::get_num_operator_effects(int op_index, bool is_axiom) const {
    (void)is_axiom;
	return op_eff[op_index].size();
}

FactPair PiMCompiledTask::get_operator_effect(int op_index, int eff_index, bool is_axiom) const {
    (void)is_axiom;
    return op_eff[op_index][eff_index];
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
    std::vector<int> new_values(domain_size.size(), 0);

    for (const auto &[fact_pairs, index] : meta_atom_map) {
        const FactPair &first_atom = fact_pairs.first;
        const FactPair &second_atom = fact_pairs.second;

        bool first_valid = (first_atom.var == -1 || values[first_atom.var] == first_atom.value);
        bool second_valid = (second_atom.var == -1 || values[second_atom.var] == second_atom.value);

        if (first_valid && second_valid) {
            new_values[index] = 1;
        }
    }

    values = std::move(new_values);
}

std::shared_ptr<AbstractTask> build_pi_m_compiled_task(
	const shared_ptr<AbstractTask> &parent) {
	return make_shared<PiMCompiledTask>(parent);
}
}