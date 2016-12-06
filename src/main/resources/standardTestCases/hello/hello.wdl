task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
    sleep 2
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
