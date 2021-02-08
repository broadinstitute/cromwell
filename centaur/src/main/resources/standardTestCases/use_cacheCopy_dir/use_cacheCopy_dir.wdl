version 1.0

task hello {

  command {
    echo "Hello again" > hello.txt
  }
  output {
    File hello_file = "hello.txt"
  }
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
}

workflow wf_hello {
  call hello
  output {
     File out = hello.hello_file
  }
}
