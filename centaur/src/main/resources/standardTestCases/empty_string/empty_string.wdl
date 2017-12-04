task hello {
  command {
    echo "Hello!"
  }
  output {
    String empty = ""
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task goodbye {
  String emptyInputString
  command {
    echo "${emptyInputString}"
  }
  output {
    String empty = read_string(stdout())
  }
  runtime {
   docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  call hello
  call goodbye {input: emptyInputString=hello.empty }
  output {
   hello.empty
   goodbye.empty
  }
}
