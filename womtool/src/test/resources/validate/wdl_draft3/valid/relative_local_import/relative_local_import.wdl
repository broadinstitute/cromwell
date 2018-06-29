version 1.0

import "../task_only/task_only.wdl" as imports

workflow relative_local_import {

  input {
    File in_file
  }

  call imports.task_only {
    input:
      in_file = in_file
  }
}
