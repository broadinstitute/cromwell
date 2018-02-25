task notroot {
  command {
    whoami
  }

  runtime {
    docker: "mcovarr/notroot:v1"
  }

  output {
    String user = read_string(stdout())
  }
}

workflow woot {
  call notroot

  output {
    String notrootUser = notroot.user
  }
}
