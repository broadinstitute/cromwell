version 1.0

# A chain of out-of-order declarations:

workflow foo {

  Int y = x
  Int a = 5
  Int x = a

  output {
    Int z = y
  }
}
