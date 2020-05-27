version 1.0

workflow multiple_overridable_optionals {
  call foo
  call foo as foo_a

  meta {
    allowNestedInputs: true
  }
}

task foo {
  input {
    Int? x
    Int y = 5
  }
  command { }
  runtime { docker: "ubuntu:latest" }
  output {
    Int? z = x
  }
}
