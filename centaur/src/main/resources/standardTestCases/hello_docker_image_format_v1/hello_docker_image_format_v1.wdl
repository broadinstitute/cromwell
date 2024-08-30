task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:utopic-20150319"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
