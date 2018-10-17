version 1.0

workflow bad_assignment_type {
  input {
    Array[String] array_of_strings
  }

  call needs_an_int { input: i = array_of_strings }

  output {
    Int j = needs_an_int.j
  }
}

task needs_an_int {
  input {
    Int i
  }
  command { }
  output {
    Int j = i
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
