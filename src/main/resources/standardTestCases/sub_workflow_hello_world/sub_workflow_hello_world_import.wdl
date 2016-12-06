task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
    sleep 2
  }
  runtime {
      docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string(stdout())
  }
}

workflow wf_hello {
  String wf_hello_input = "world"
  
  call hello {input: addressee = wf_hello_input }
  
  output {
    String salutation = hello.salutation
  }
}