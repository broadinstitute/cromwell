version draft-3

workflow needs_an_input {
  meta {
    description: "Equivalent to the invalid 'missing input' case, except that we don't provide an inputs JSON to validate against."
  }
  input {
    Int x
  }
}
