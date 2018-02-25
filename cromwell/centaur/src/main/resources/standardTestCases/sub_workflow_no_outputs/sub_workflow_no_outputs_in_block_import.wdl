
workflow subhello {
  Array[String] greeting_pieces

  call hello {
    input: inputs = greeting_pieces
  }

  String salutation_length = length(hello.out)

  # Neither the call output nor the declaration should be considered outputs to the subworkflow
  output { }
}

task hello {
  Array[String] inputs

  command {}

  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Array[String] out = inputs
  }
}
