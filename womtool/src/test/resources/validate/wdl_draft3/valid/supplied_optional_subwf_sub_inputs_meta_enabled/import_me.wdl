version 1.0

workflow sub_wf {
  input {
    Int y
  }
  # Calls foo but doesn't provide an 'x'. That's fine because the input was optional, but now the outer WF cannot override it.
  call foo { input: y = y }

  # No outputs, but we don't want that to be the error:
  output { }
  meta {allowNestedInputs: true}
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
