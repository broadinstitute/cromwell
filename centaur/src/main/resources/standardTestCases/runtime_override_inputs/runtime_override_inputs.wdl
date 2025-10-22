version development-1.1

task hello {
  input {
    String addressee
  }

  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
    maxRetries: 3
    preemptible: 1
  }
}

workflow wf_hello {
  call hello
  output {
     String salutation = hello.salutation
  }
}
