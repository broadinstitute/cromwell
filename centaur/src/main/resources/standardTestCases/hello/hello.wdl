task hello {
  command {
    echo "Hello DRAGEN!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "dragenpublicwestus.azurecr.io/dragen-xrt/el7/x86_64:3.7.8m"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
