version 1.0

import "subworkflow.wdl" as sub

workflow subworkflow_input {
  call sub.subwf
}
