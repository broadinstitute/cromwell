import "sub_workflow_hello_world.wdl" as hello
import "forkjoin.wdl"
import "sub_workflow_interactions.wdl" as interactions

task join {
  Int grepCount
  Int wcCount
  command {
    expr ${wcCount} / ${grepCount}
  }
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow recursive_imports {
  call hello.sub.wf_hello {}
  call forkjoin.join { input: wcCount=15, grepCount=10 }
  call interactions.counter.countTo { input: value = 3 }

  output {
    wf_hello.*
    join.*
  }
}
