version development-1.1

import "sub_runtime_override_inputs.wdl" as sub_runtime

task hello {
  input {
    String addressee
  }

  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
    maxRetries: 3
    preemptible: 1
  }
}

workflow wf_hello {

  meta {
    allowNestedInputs: true
  }

  call hello
  call sub_runtime.wf as sub_runtime_call

  output {
     String salutation = hello.salutation
  }
}
