task task_only {

  File in_file

  command {
    wc -w < ${in_file}
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    String word_count = read_int(stdout())
  }
}
