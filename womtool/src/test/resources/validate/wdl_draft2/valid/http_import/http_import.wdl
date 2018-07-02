import "https://raw.githubusercontent.com/broadinstitute/cromwell/5e0197d1c016d4c802ef3c2890f0ca4e0ca542c1/womtool/src/test/resources/validate/wdl_draft2/valid/task_only/task_only.wdl" as imports

workflow http_import {

  File in_file

  call imports.task_only {
    input:
      in_file = in_file
  }
}
