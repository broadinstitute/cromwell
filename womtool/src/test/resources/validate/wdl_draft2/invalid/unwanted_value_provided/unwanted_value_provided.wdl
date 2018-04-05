workflow unwanted_value_provided {
  call has_unoverridable_value { input: x = 5, y = 6 }
}

task has_unoverridable_value {
  # x is suppliable:
  Int x
  # y is not (it has an upstream dependency):
  Int y = x

  command { echo ${y} }

  output {
    Int z = y
  }
}
