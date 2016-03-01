task failOnStderr {
    command <<<
        echo "OH NO!" >&2
    >>>
    output {
        String ohno = read_string(stderr())
    }
    runtime {
        docker: "ubuntu:latest"
        failOnStderr: true
    }
}

workflow runtime_failOnStderr {
    call failOnStderr
}
