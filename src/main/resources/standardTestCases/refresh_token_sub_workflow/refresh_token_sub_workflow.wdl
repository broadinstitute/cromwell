import "<<IMPORT>>" as sub

task hello {
  File iFile
  String addressee = read_string(iFile)
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow refresh_token_sub_workflow {
  call hello
  call sub.wf_hello
  output {
     hello.salutation
     wf_hello.salutation
  }
}
