task dockerhub {
    command {
        echo "bonjour tout le monde !"
    }
    runtime {
        docker: "ubuntu:precise-20161209"
        backend: "LocalNoDocker"
    }
}

workflow backend_with_no_docker {
    call dockerhub
}
