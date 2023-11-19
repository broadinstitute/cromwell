version 1.0

import "testworkflow.wdl" as sub

task start {
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
  command {
    echo "${2/2}"
  }
  output {
    String out = read_string(stdout())
  }
}

task done {
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
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
