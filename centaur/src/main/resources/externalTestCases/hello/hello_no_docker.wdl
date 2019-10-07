task hello {
  String addressee = "you"
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
}

workflow wf_hello {
  call hello
  output {
     hello.salutation
  }
}
