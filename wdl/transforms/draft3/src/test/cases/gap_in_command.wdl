version 1.0

workflow my_workflow {
  call my_task
}

task my_task {
  command {
    echo "hi"

    echo "bye"
  }

  output {
    Array[String] lines = read_lines(stdout())
  }
}
