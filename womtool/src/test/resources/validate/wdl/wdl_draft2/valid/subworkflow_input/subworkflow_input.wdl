import "subworkflow.wdl" as sub

workflow subworkflow_input {
  call sub.subwf
}
