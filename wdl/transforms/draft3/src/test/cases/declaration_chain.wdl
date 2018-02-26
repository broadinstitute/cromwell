version draft-3

# A chain of out-of-order declarations:

workflow foo {

  Int y = x
  Int a = 5
  Int x = a

  output {
    Int z = y
  }
}
