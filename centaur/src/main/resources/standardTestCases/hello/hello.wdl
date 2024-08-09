task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
  meta {
    volatile: "true"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
