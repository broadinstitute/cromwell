version draft-3

workflow unwanted_value_provided {
  call has_unoverridable_value { input: x = 5, y = 6 }
}

task has_unoverridable_value {
  input {
    Int x
  }

  Int y = x

  command { echo ~{y} }

  output {
    Int z = y
  }
}
