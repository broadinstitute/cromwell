task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
}

workflow local_backend {
  call hello
  output {
     hello.salutation
  }
}
