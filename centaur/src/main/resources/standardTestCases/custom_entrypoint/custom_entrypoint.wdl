version 1.0

workflow custom_entrypoint_wf {
    call custom_entrypoint_task
    output {
        String hello = custom_entrypoint_task.hello
    }
}

task custom_entrypoint_task {
    meta {
        # Don't allow call-caching
        volatile: true
    }
    command <<<
        echo hello world
    >>>
    runtime {
        # This image has a pre-defined entrypoint so we _have_ to override it:
        docker: "broadinstitute/cromwell-docker-test:custom_entrypoints"
    }
    output {
        String hello = read_string(stdout())
    }
}
