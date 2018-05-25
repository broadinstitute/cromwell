workflow unexpected_call_input_failure {

  call hello {input:
    addressee = "mon ami", # This is an expected input
    greeting = "bonjour"   # But this is an unexpected input! Should be an error!
  }
  
  output {
    String salutation = hello.salutation
  }
}

task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string(stdout())
  }
}
