task quay {
    command {
        echo "hello"
    }
    runtime {
        docker: "quay.io/fedora/fedora-minimal@sha256:63f785bf6185a63332fd9ccdaaabc66f47bc63564e2cecc6f99194002ca73ac3"
    }
}

workflow docker_hash_quay {
    call quay 
}
