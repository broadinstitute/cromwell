version 1.0

workflow foo {
  input {
    Int x
  }

  Int y = x

  output {
    Int z = y
  }
}
