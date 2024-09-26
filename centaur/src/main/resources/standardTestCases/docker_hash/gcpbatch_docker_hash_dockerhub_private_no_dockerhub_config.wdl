task dockerhub {
    command {
        echo "hello"
    }
    runtime {
        docker: "broadinstitute/cloud-cromwell:dev"
        backend: "GCPBATCHNoDockerHubConfig"
    }
}

workflow docker_hash_dockerhub_private {
    call dockerhub
}
