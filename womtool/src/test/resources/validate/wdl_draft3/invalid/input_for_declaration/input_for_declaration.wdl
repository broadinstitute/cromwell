version draft-3

workflow input_for_declaration {
  meta {
    description: "This submission is invalid because it tries to input a value which is an intermediate declaration"
  }

  input {

  }

  # Not an input:
  Int x = 5
}
