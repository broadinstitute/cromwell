task dockerhub {
    command {
        echo "hello"
    }
    runtime {
        docker: "broadinstitute/cloud-cromwell:2024-08-31"
        backend: "GCPBATCHGoogleSecretManager"
    }
}

workflow docker_hash_dockerhub_private {
    call dockerhub
}
