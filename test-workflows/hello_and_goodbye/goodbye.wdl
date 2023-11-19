version 1.0

task sayFarewell {
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