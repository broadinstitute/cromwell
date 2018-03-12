task notroot {
  command {
    echo $HOME
  }

  runtime {
    docker: "broadinstitute/cromwell-docker-test:notroot"
  }

  output {
    String home = read_string(stdout())
  }
}

workflow woot {
  call notroot

  output {
    String notrootHome = notroot.home
  }
}
