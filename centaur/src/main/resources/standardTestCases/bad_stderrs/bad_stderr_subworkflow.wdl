version 1.0

import "bad_stderr_mid_workflow.wdl"

workflow calls_bad_stderr_subworkflow {
    scatter(i in range(1000)) {
        call bad_stderr_mid_workflow.bad_stderr
    }
}
