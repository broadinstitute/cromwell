task gcr {
    command {
        echo "hello"
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-cromwell-dev/centaur-ubuntu:v0"
    }
}

workflow docker_hash_gcr_private {
    call gcr
}
