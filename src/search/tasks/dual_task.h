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
};


std::shared_ptr<AbstractTask> build_dual_task(
	const std::shared_ptr<AbstractTask> &parent);

}

#endif
