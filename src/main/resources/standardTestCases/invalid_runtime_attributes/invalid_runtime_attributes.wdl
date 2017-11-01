task invalid_runtime_attributes {
    command { # NOOP }
    runtime {
        docker: "ubuntu:latest"
        continueOnReturnCode: "oops"
    }
}

workflow invalid_runtime_attributes_wf {
    call invalid_runtime_attributes
}
