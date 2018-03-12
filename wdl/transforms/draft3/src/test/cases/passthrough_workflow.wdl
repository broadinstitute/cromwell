version draft-3

workflow foo {
  input {
    Int x
  }

  Int y = x

  output {
    Int z = y
  }
}
