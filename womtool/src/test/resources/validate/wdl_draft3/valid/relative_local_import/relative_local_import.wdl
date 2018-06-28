version 1.0

import "./womtool/src/test/resources/validate/wdl_draft3/valid/task_only/task_only.wdl" as imports

workflow http_import {

  input {
    File in_file
  }

  call imports.task_only {
    input:
      in_file = in_file
  }
}