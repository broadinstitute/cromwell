task gcr {
    command {
        echo "hello"
    }
    runtime {
        docker: "gcr.io/google-containers/ubuntu-slim:0.14"
    }
}

workflow docker_hash_gcr {
    call gcr
}
