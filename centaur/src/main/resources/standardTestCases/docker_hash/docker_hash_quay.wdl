task quay {
    command {
        echo "hello"
    }
    runtime {
        docker: "quay.io/broadinstitute/cromwell-tools:v2.4.1"
    }
}

workflow docker_hash_quay {
    call quay 
}
