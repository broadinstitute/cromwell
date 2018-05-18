version 1.0

workflow bad_assignment_type {
  input {
    Array[String] y = 57
  }
  output {
    Array[String] strings = y
  }
}
