import "sub_workflow_hello_world_import.wdl" as sub

workflow main_workflow {
    call sub.wf_hello { input: wf_hello_input = "sub world" }
    output {
        String main_output = wf_hello.salutation
    }
}