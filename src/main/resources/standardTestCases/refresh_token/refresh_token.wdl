task hello {
  File iFile
  String addressee = read_string(iFile)
  command {
    echo "Hello ${addressee}!"
    sleep 2
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow refresh_token {
  call hello
  output {
     hello.salutation
  }
}
