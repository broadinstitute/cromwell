task mkdir {
    command {
        mkdir empty_dir
    }
    output {
        File empty_dir = "empty_dir"
    }
    runtime {
        docker: "manifoldai/alpine-plus-bash:3.20.0"
    }
}

workflow docker_alpine {
    call mkdir
}
