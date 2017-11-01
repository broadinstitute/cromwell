task dockerhub {
    command {
        echo "hello"
    }
    runtime {
        docker: "broadinstitute/cloud-cromwell:dev"
    }
}

workflow docker_hash_dockerhub_private {
    call dockerhub
}
