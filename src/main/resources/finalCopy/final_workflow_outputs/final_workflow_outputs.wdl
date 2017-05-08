task hello {
  command {
    echo "Hello !" > test.out
  }
  output {
    File out = "test.out"
  }
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.out
  }
}
