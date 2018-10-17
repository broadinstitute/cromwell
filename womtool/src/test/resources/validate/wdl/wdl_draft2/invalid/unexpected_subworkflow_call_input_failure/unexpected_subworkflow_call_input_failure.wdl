import "unexpected_subworkflow_call_input_failure_import.wdl" as subwf

workflow unexpected_subworkflow_call_input_failure {
  call subwf.subwf { input: i = 10, j = 20 }
}
