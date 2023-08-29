task hello {
  command {
    echo "Hello DRAGEN!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "whatever38fj.azurecr.io/f00:v6.456"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
