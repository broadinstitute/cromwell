version 1.0

workflow overridable_optionals {
  call foo as foo0

  scatter(i in []) {
    call foo as foo1

    scatter(j in []) {
      call foo as foo2
    }
  }
}

task foo {
  input {
    Int z = 5
  }
  command { }
  runtime { docker: "ubuntu:latest" }
  output {
    Int? out = z
  }
}
