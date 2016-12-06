task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
    sleep 2
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
   docker: "/fauxbuntu:nosuchversion"
  }
}

workflow wf_hello {
  call hello
}
