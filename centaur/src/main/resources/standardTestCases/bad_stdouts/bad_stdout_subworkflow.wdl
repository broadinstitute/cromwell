version 1.0

import "bad_stdout_mid_workflow.wdl"

workflow calls_bad_stdout_subworkflow {
    scatter(i in range(1000)) {
        call bad_stdout_mid_workflow.bad_stdout
    }
}
