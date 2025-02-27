#ifndef DUAL_TASK_H
#define DUAL_TASK_H

#include "delegating_task.h"

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>

class AbstractTask;

namespace extra_tasks {
class DualTask : public tasks::DelegatingTask {
    struct FactPairHash {
        size_t operator()(const FactPair &fact) const {
            return fact.var * 100003 + fact.value;
        }
    };
    std::map<FactPair, int> strips_atom_map;
  	std::vector<int> initial_state_values;
	std::vector<FactPair> goals;
	std::vector<std::vector<FactPair>> dual_operator_pre;
    std::vector<std::vector<FactPair>> dual_operator_eff;

    public:
		DualTask(const std::shared_ptr<AbstractTask> &parent);

    void init_strips_variables();
    void setup_init_and_goal_states();
    void setup_dual_operators();

    void print_task();
    const std::vector<FactPair> get_goals() const;

    virtual int get_num_variables() const override;
    virtual int get_variable_domain_size(int var) const override;
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
    virtual void convert_state_values_from_parent(
        std::vector<int> &values) const override;

};


std::shared_ptr<AbstractTask> build_dual_task(
	const std::shared_ptr<AbstractTask> &parent);

}

#endif
