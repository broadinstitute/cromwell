task dockerhub {
    command {
        echo "hello"
    }
    runtime {
        docker: "ubuntu:precise-20161209"
    }
}

workflow docker_hash_dockerhub {
    call dockerhub
}
