version 1.0

task sayFarewell {
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
  input {
    String farewell
  }
  command {
    echo "${farewell}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow goodbye {
  input {
    String farewell
  }
  call sayFarewell {
    input: farewell = farewell
  }
  output {
    String out = sayFarewell.out
  }
}