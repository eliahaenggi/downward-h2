#ifndef PI_M_COMPILED_TASK_H
#define PI_M_COMPILED_TASK_H

#include "delegating_task.h"

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class AbstractTask;

namespace extra_tasks {
class PiMCompiledTask : public tasks::DelegatingTask {

    struct FactPairHash {
        size_t operator()(const FactPair &fact) const {
            return fact.var * 100003 + fact.value;
        }
    };
	struct AtomPairHash {
		std::size_t operator()(std::pair<FactPair, FactPair> pair) const {
            const int MOD = 100003;
            std::size_t h1 = pair.first.var * MOD + pair.first.value;
            std::size_t h2 = pair.second.var * MOD + pair.second.value;
            return h1 * MOD + h2;
        }
	};

    std::unordered_map<std::pair<FactPair,FactPair>, int, AtomPairHash> meta_atom_map;
	std::vector<int> initial_state_values;
	std::vector<FactPair> goals;


    struct MetaOperator {
    int parent_id;
    FactPair s_atom;
    std::vector<FactPair> preconditions;
    std::vector<FactPair> effects;
    int cost;

    MetaOperator(const int par_id, const FactPair s,
                 const std::vector<FactPair> pre, const std::vector<FactPair> eff, const int c) :
				parent_id(par_id), s_atom(s), preconditions(pre), effects(eff), cost(c) {}
	};

    std::vector<MetaOperator> meta_operators;
public:
	PiMCompiledTask(const std::shared_ptr<AbstractTask> &parent);
    // Main functions to set up compiled task
    void init_meta_atom_map();
    void setup_init_and_goal_states();
    void setup_meta_operators();

    // Helper functions
    FactPair translate_into_meta_atom(FactPair first_atom, FactPair second_atom);
    bool contradict_precondition(int op_id, FactPair s_atom);
    std::vector<FactPair> generate_meta_preconditions(int op_id);
    std::vector<FactPair> generate_meta_effects(int op_id, std::unordered_set<int> &effect_vars);

    // Functions to access compiled task transformation
    virtual int get_num_variables() const override;
    virtual int get_variable_domain_size(int var) const override;
    virtual int get_operator_cost(int index, bool is_axiom) const override;
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


std::shared_ptr<AbstractTask> build_pi_m_compiled_task(
	const std::shared_ptr<AbstractTask> &parent);

}

#endif //PI_M_COMPILED_TASK_H
