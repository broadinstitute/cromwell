version 1.0

workflow passthrough_value {
  input {
    Int x
  }

  Int y = x

  output {
    Int z = y
  }
}
