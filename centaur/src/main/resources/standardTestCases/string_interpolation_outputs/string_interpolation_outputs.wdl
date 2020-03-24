version 1.0

task what_is_one {
  Int one = 1

  command {
    echo "one" > 1.txt
  }
  output {
    String one_string = "one is ~{read_string("1.txt")}"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow string_interpolation_outputs {
  call what_is_one

  output {
    String one_is = what_is_one.one_string
  }
}

