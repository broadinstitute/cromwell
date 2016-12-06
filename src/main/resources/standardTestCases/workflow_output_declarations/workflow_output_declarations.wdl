task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
    sleep 2
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string(stdout())
  }
}

workflow wf_output_declarations {
  Array[String] arr = ["Me"]
  call hello {input: addressee = "You" }
  
  scatter(i in arr) {
    call hello as hello2 {input: addressee = i }
  }
  
  output {
    String out1 = hello.salutation
    Array[String] out2 = hello2.salutation
    String out3 = out2[0] + out1
  }
}