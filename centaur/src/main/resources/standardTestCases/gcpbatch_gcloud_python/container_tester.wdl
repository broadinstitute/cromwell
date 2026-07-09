version 1.0

workflow container_tester {
    input {
        String container_image
    }

    call run_container {
        input:
            container_image = container_image
    }
}

task run_container {
    input {
        String container_image
    }

    # `gcloud` should successfully find Python, run, and exit (AN-601)
    command <<<
        printenv
        gcloud --version
    >>>

    runtime {
        docker: container_image
    }

    meta {
        volatile: true
    }
}
