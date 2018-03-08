version draft-3

workflow simple_conditional {
  Boolean bool = true

  if (bool) {
    Int j = i + 10
  }

  output {
    Int? j_maybe = j
  }
}
