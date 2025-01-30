#include "pi_m_compiled_task.h"

#include "../task_utils/task_properties.h"



using namespace std;
namespace extra_tasks {

PiMCompiledTask::PiMCompiledTask(const shared_ptr<AbstractTask> &parent) : DelegatingTask(parent) {
  	TaskProxy parent_proxy(*parent);
	if (task_properties::has_conditional_effects(parent_proxy)) {
		ABORT("DomainAbstractedTask doesn't support conditional effects.");
	}
}

shared_ptr<AbstractTask> PiMCompiledTask::get_task() const {
	return task;
}

std::shared_ptr<AbstractTask> build_pi_m_compiled_task(
	const shared_ptr<AbstractTask> &parent) {
	return PiMCompiledTask(parent).get_task();
}
}