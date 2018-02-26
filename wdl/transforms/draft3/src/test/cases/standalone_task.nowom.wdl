version draft-3

task standalone {
  input {
    String bar
  }

  output {
    String out = bar
  }

  command {
    echo ${bar}
  }

  runtime {
    docker: "someFakeDockerRuntime"
  }
}
