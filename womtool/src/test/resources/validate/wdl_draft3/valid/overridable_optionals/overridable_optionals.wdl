version 1.0

workflow overridable_optionals {
  input {
    Int x
  }
  call foo {input: x=x}
  meta {
    allowNestedInputs: true
  }
}

task foo {
  input {
    Int x
    Int? y
    Int z = 5
  }
  command { }
  runtime { docker: "ubuntu:latest" }
  output {
    Int? out = x
  }
}
