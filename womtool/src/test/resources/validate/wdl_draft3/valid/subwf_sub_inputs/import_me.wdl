version 1.0

workflow sub_wf {
  input {
    Int y
  }
  # Calls foo but doesn't provide an 'x'. That's fine because the inputs are optional
  call foo { input: y = y }

  output {
    Int z = foo.z
  }
}

task foo {
  input {
    Int? x
    Int y
  }
  command {
  }
  output {
    Int z = y
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
