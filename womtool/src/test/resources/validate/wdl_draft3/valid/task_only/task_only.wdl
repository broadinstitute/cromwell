version 1.0

task task_only {
  input {
    File in_file
  }
  command {
    wc -w < ~{in_file}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String word_count = read_int(stdout())
  }
}
