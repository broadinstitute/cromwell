version 1.0

workflow placeholder_attributes {
  call uses_attributes
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
    echo ~{tre="good" false="bad" t}
    echo ~{fale="bad" true="good" f}
    echo ~{defult="bad" supplied}
    echo ~{seep="___" items}
    echo ~{true="good" t}
    echo ~{false="good" f}
    echo ~{true="good" true="what" false="bad" t}
    echo ~{false="bad" false="what" true="good" f}
    echo ~{default="bad" default="even worse" supplied}
    echo ~{sep="___" sep="uuuuh" items}
  }

  runtime {
    docker: "ubuntu:latest"
  }
}

