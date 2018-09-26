version 1.0

import "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/centaur/src/main/resources/standardTestCases/wdl_draft3/draft3_infer_version/draft3_infer_version.wdl" as infered

workflow draft3_http_imports {
  call infered.draft3_infer_version

  output {
    Int j = draft3_infer_version.j
  }
}
