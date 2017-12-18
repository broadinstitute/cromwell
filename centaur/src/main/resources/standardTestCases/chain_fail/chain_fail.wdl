import "chain_fail_import.wdl" as sub

workflow chain_fail {
  call fail
    
  Array[Boolean] is = [ fail.o ]

  # This scatter will never be runnable because the fail task will... fail
  scatter (j in is) {
    call fail as scatter_fail
  }
  
  if (fail.o) {
    call fail as conditional_fail
  }
  
  call sub.wf_hello { input: wf_hello_input = fail.o }

  output {
    Array[Boolean] scatter_o = scatter_fail.o
    Boolean? conditional_o = conditional_fail.o
    String sub_workflow_o = wf_hello.salutation
  }
}


task fail {
  command {
    echo false
    exit 1
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Boolean o = read_boolean(stdout())
  }
}
