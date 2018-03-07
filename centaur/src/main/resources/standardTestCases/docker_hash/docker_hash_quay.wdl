task quay {
    command {
        echo "hello"
    }
    runtime {
        docker: "quay.io/broadinstitute/cromwell-docker-test:centaur"
    }
}

workflow docker_hash_quay {
    call quay 
}
