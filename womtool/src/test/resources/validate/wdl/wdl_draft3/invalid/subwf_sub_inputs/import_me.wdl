version 1.0

workflow sub_wf {
  # Calls foo but doesn't provide an 'x', so this workflow shouldn't (strictly) be allowed as a sub-workflow
  call foo
}

task foo {
  input {
    Int x
  }
  command {
  }
  output {
    Int y = x
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
