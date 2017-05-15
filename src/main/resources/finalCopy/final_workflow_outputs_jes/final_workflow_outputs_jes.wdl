task hello {
  command {
    echo "Hello !" > test.out
  }
  output {
    File out = "test.out"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.out
  }
}
