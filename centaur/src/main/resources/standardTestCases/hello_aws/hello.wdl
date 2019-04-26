task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "952500931424.dkr.ecr.us-east-1.amazonaws.com/broadinstitute/cromwell:latest"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
