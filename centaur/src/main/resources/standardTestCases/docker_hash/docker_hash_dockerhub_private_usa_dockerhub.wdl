task dockerhub {
    command {
        echo "hello"
    }
    runtime {
        docker: "broadinstitute/cloud-cromwell:dev"
        backend: "Papiv2USADockerhub"
    }
}

workflow docker_hash_dockerhub_private {
    call dockerhub
}
