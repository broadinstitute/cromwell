version 1.0

import "import_me.wdl" as imported

workflow bad_nested_task_output_assignment_type {

  if (true) {
    scatter(x in range(4)) {
      call imported.inner_workflow
    }
  }

  output {
    String bad_gather = inner_workflow.x
  }
}
