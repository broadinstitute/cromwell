version 1.0



workflow simple_conditional {
  Boolean bool = true

  Int i = 5

  if (bool) {
    Int j = i + 10
  }

  output {
    Int? j_maybe = j
  }
}
