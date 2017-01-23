task nobody {
  command {
    echo $HOME
  }

  runtime {
    docker: "ubuntu:latest"
    docker_user: "nobody"
  }

  output {
    String home = read_string(stdout())
  }
}

workflow woot {
  call nobody

  output {
    String nobodyHome = nobody.home
  }
}
