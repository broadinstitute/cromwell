task prefix {
  command {
    echo hello \
    | wc -l
  }
  output {
    String out = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_prefix {
  call prefix
  output {
     String prefix_out = prefix.out
  }
}
