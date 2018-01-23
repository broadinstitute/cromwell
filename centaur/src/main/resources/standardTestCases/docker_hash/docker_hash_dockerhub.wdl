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

    Array[String] hostnames = ["", "docker.io/", "index.docker.io/"]

    scatter (hostname in hostnames) {
        call dockerhub { input: hostname = hostname }
    }
}

