task hello {
  command {
    echo "Hello " > test.out
  }
  output {
    File out = "test.out"
  }
}

workflow wf_hello {
  call hello
  output {
     hello.out
  }
}
