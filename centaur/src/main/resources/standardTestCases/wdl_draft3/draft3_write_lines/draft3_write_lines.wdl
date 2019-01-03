version 1.0

# Testing write_lines() within a command block
workflow draft3_write_lines {
  call a2f {
    input: strings = ["a","b","c","d"]
  }

  call f2a {
    input: i=a2f.out
  }

  call a2f as a2f_second {
    input: strings = f2a.out
  }

  output {
    String x_out = a2f_second.x
  }
}

task f2a {
  input {
    File i
  }

  command {
    cat ~{i}
  }

  output {
    Array[String] out = read_lines(stdout())
  }

  runtime {
    docker:"ubuntu:latest"
  }
}

task a2f {
  input {
    Array[String] strings
  }

  command {
    cat ~{write_lines(strings)}
  }

  output {
    File out = stdout()
    String x = read_string(out)
  }

  runtime {
    docker:"ubuntu:latest"
  }
}
