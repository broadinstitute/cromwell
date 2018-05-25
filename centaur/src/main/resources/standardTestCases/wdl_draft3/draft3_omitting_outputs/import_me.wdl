version 1.0

task foo {
  command {
    echo hi
  }
  output {
    String x = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow foo_in_a_box {
  call foo
  output {
    String x = foo.x
  }
}
