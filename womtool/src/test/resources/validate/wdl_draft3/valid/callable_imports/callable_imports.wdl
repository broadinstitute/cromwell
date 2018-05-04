version 1.0

import "my_task.wdl"
import "my_workflow.wdl" as mai_workflau

workflow callable_imports {
  call my_task.my_task {
    input: s = "foo"
  }

  call mai_workflau.my_workflow as may_wurkfloo {
    input: i = my_task.i
  }

  output {
    Int i = my_task.i
    Int o = may_wurkfloo.o
  }
}
