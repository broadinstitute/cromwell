version 1.0

workflow Test {

  Int a = 5

  # We are checking that b is of type Int? and not Int?? or Int??? due to nesting
  if (true) {
    if (true) {
      if (true) {
        Int b = 5
      }
    }
  }

  Int c = select_first([a, b])
}
