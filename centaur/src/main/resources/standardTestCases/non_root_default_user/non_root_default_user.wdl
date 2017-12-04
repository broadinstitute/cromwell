task notroot {
  command {
    echo $HOME
  }

  runtime {
    docker: "mcovarr/notroot:v1"
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
