#ifndef PI_M_COMPILED_TASK_H
#define PI_M_COMPILED_TASK_H

#include "delegating_task.h"

#include <memory>
#include <string>
#include <vector>

class AbstractTask;

namespace extra_tasks {
class PiMCompiledTask : public tasks::DelegatingTask {
	std::vector<int> domain_size;
	std::vector<int> initial_state_values;
	std::vector<FactPair> goals;
	std::vector<std::vector<std::string>> fact_names;
	std::vector<std::vector<int>> value_map;
	std::shared_ptr<AbstractTask> task;
public:
	PiMCompiledTask(const std::shared_ptr<AbstractTask> &parent);

    std::shared_ptr<AbstractTask> get_task() const;

};


std::shared_ptr<AbstractTask> build_pi_m_compiled_task(
	const std::shared_ptr<AbstractTask> &parent);

}

#endif //PI_M_COMPILED_TASK_H
