task mkdir {
    command {
        mkdir empty_dir
    }
    output {
        File empty_dir = "empty_dir"
    }
    runtime {
        docker: "python:alpine"
    }
}

workflow docker_alpine {
    call mkdir
}
