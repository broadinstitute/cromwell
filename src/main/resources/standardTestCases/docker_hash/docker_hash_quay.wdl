task quay {
    command {
        echo "hello"
    }
    runtime {
        docker: "quay.io/aptible/ubuntu:12.04"
    }
}

workflow docker_hash_quay {
    call quay 
}
