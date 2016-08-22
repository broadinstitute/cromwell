task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
   docker: "/fauxbuntu:nosuchversion"
  }
}

workflow hello {
  call hello
}
