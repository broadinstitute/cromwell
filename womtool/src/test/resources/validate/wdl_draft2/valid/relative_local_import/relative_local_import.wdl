import "./womtool/src/test/resources/validate/wdl_draft2/valid/task_only/task_only.wdl" as imports

workflow http_import {

  File in_file

  call imports.task_only {
    input:
      in_file = in_file
  }
}