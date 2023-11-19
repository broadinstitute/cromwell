version 1.0

import "speak.wdl" as sub

task start {
  command {
    echo "Starting...."
  }
  output {
    String out = read_string(stdout())
  }
}

task done {
  command {
    echo "Done"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow main_workflow {
  call start
  call sub.speak
  call done
}