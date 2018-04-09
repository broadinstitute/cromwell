version draft-3

workflow missing_input {
  meta {
    description: "Equivalent to the valid 'needs_an_input' case, except that we provide an inputs JSON which is insufficient."
  }
  input {
    Int x
  }
}
