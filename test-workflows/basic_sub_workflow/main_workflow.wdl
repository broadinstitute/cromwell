version 1.0

import "testworkflow.wdl" as sub

task start {
  command {
    echo "${2/2}"
  }
  output {
    String out = read_string(stdout())
  }
}

task done {
  command {
    echo "${2/0}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow main_workflow {
  call start
  call sub.wf_scattering
  call done
}
