version 1.0

import "sub_workflow_delete_import.wdl" as sub_workflow_delete_import

workflow sub_workflow_delete {
    call sub_workflow_delete_import.sub_workflow_delete_import as sub_call {
        input:
    }

    output {
        File keep = sub_call.keep
    }
}
