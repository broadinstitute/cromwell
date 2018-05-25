version 1.0

workflow use_previous_file_outputs {
  call read_previous_string
  output {
    String hi_out = read_previous_string.hi_out
  }
}

task read_previous_string {
  command {
    echo hi > output.txt
  }
  output {
    File hi_out_file = "output.txt"
    String hi_out = read_string(hi_out_file)
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
