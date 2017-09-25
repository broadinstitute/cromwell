
workflow subhello {
  Array[String] greeting_pieces

  call hello {
    input: inputs = greeting_pieces
  }

  # Confirm referencing call outputs from subworkflows does not fail validation.
  String salutation_length = length(hello.out)

  output {}
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
