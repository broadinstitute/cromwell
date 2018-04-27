import "subworkflow.wdl" as sub

workflow missing_subworkflow_input {
  call sub.subwf
}
