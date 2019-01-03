task notroot {
  String empty = ""
  command {
    echo $HOME ${empty}
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
