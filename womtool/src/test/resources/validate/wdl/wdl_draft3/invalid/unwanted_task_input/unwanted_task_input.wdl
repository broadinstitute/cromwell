version 1.0

workflow unwanted_value_provided {
  call has_unoverridable_value { input: x = 5, z = 6 }
}

task has_unoverridable_value {
  input {
    Int x
    Int y = x
  }

  Int z = x

  command { echo ~{y} }

  output {
    Int z_out = z
  }
}
