version 1.0

import "https://raw.githubusercontent.com/broadinstitute/cromwell/5e0197d1c016d4c802ef3c2890f0ca4e0ca542c1/womtool/src/test/resources/validate/wdl_draft3/valid/task_only/task_only.wdl" as imports

workflow http_import_missing_input {
  call imports.task_only
}
