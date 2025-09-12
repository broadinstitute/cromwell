version 1.0

workflow standalone_task {
    call standalone { input: bar="Hello, World!" }
}

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
    container: "someOtherDockerRuntime"
  }
}
