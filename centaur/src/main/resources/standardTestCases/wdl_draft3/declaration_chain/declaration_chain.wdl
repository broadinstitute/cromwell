version 1.0

# A chain of out-of-order declarations:

workflow declaration_chain {

  Int y = x
  Int a = 5
  Int x = a

  output {
    Int z = y
  }
}
