version draft-3

workflow standalone_task {}

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
