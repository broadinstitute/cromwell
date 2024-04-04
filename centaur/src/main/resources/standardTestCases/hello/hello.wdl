task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "terrabatchdev.azurecr.io/nginx:v1"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
