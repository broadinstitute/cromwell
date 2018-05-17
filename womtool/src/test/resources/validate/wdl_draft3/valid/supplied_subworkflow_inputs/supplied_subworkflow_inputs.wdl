version 1.0

import "subworkflow.wdl" as sub

workflow subworkflow_input {
  call sub.subwf { input:
    i = 5,
    j = 10,
    k = 15
  }
}
