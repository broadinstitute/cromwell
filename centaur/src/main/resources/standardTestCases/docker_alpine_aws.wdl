task mkdir {
    command {
        mkdir empty_dir
    }
    output {
        File empty_dir = "empty_dir"
    }
    runtime {
        docker: "alpine:3.5"
    }
}

workflow docker_alpine {
    call mkdir
}
