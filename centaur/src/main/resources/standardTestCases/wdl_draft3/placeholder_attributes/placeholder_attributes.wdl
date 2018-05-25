version 1.0

workflow placeholder_attributes {
  call uses_attributes

  output {
    Array[String] lines = uses_attributes.lines
  }
}

task uses_attributes {

  Boolean t = true
  Boolean f = true
  Array[Int] items = [1, 2, 3, 4, 5]

  String? supplied = "good"

  input {
    String? unsupplied
  }

  command {
    echo ~{true="good" false="bad" t}
    echo ~{false="bad" true="good" f}
    echo ~{default="bad" supplied}
    echo ~{default="good" unsupplied}
    echo ~{sep="___" items}
  }

  output {
    Array[String] lines = read_lines(stdout())
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
