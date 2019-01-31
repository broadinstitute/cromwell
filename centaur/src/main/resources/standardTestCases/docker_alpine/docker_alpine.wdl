task mkdir {
    command {
        mkdir empty_dir
    }
    runtime {
        docker: "alpine:3.5"
        backend: "LocalBourneShell"
    }
}

workflow docker_alpine {
    call mkdir
}
