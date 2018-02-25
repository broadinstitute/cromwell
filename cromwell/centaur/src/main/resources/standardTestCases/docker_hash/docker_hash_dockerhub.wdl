task dockerhub {
    String hostname
    command {
        echo "hello"
    }
    runtime {
        docker: hostname + "ubuntu:precise-20161209"
    }
}

workflow docker_hash_dockerhub {

    call dockerhub as dockerWithoutHost { input: hostname = "" }

    call dockerhub as implicitRegistryHost { input: hostname = "docker.io/" }

    call dockerhub as defaultRegistryHost { input: hostname = "index.docker.io/" }
}

