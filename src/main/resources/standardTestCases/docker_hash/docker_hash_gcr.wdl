task gcr {
    command {
        echo "hello"
    }
    runtime {
        docker: "gcr.io/google-containers/ubuntu:14.04"
    }
}

workflow docker_hash_gcr {
    call gcr
}
