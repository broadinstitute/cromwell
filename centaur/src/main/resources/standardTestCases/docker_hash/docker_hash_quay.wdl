task quay {
    command {
        echo "hello"
    }
    runtime {
        docker: "quay.io/fedora/fedora-minimal@sha256:8236b62386e5a4c34b6363bcf64b20147771b4e4b9c0a24f71a4e8fa5b8703f9"
    }
}

workflow docker_hash_quay {
    call quay 
}
