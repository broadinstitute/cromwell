version 1.0

workflow inner_workflow {
  call inner_task

  output {
    String x = inner_task.x
  }
}

task inner_task {
  command {
  }
  output {
    String x = "x"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
