
workflow subhello {
  Array[String] greeting_pieces

  call hello {
    input: inputs = greeting_pieces
  }

  String salutation_length = length(hello.out)

  output {
    Int sal_len = salutation_length
  }
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
