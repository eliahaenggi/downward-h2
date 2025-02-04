#ifndef PI_M_COMPILED_TASK_H
#define PI_M_COMPILED_TASK_H

#include "delegating_task.h"

#include <memory>
#include <string>
#include <vector>
#include <map>

class AbstractTask;

namespace extra_tasks {
class PiMCompiledTask : public tasks::DelegatingTask {
	std::vector<int> domain_size;
    std::map<std::pair<FactPair,FactPair>, int> meta_atom_map;
	std::vector<int> initial_state_values;
	std::vector<FactPair> goals;
	std::vector<std::vector<std::string>> fact_names;
    std::vector<std::vector<FactPair>> old_pre;
    std::vector<std::vector<FactPair>> old_eff;
    std::vector<std::vector<FactPair>> op_pre;
    std::vector<std::vector<FactPair>> op_eff;
    std::vector<int> op_cost;
    std::vector<std::pair<int, FactPair>> op_list;
	std::shared_ptr<AbstractTask> task;
public:
	PiMCompiledTask(const std::shared_ptr<AbstractTask> &parent);

    void store_old_ops();
    void init_meta_atom_map();
    void setup_init_and_goal_states();
    void setup_new_ops();

    FactPair translate_into_meta_atom(FactPair first_atom, FactPair second_atom);
    bool contradict_precondition(int op_id, FactPair s_atom);

    void dump_compiled_task();

    virtual int get_num_variables() const override;
    virtual std::string get_variable_name(int var) const override;
    virtual int get_variable_domain_size(int var) const override;
    virtual int get_operator_cost(int index, bool is_axiom) const override;
    virtual std::string get_fact_name(const FactPair &fact) const override;
    virtual std::string get_operator_name(int index, bool is_axiom) const;
    virtual int get_num_operators() const override;
    virtual int get_num_operator_preconditions(int index, bool is_axiom) const override;
	virtual FactPair get_operator_precondition(
        int op_index, int fact_index, bool is_axiom) const override;
    virtual int get_num_operator_effects(int op_index, bool is_axiom) const override;
    virtual FactPair get_operator_effect(
        int op_index, int eff_index, bool is_axiom) const override;
    virtual int get_num_goals() const override;
    virtual FactPair get_goal_fact(int index) const override;
    virtual std::vector<int> get_initial_state_values() const override;
    virtual int get_num_operator_effect_conditions(
        int op_index, int eff_index, bool is_axiom) const override;

    struct FactPairHash {
        size_t operator()(const FactPair &fact) const {
            return fact.var * 100003 + fact.value;
        }
    };
};


std::shared_ptr<AbstractTask> build_pi_m_compiled_task(
	const std::shared_ptr<AbstractTask> &parent);

}

#endif //PI_M_COMPILED_TASK_H
