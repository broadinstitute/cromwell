version 1.0

import "subworkflow.wdl" as sub

workflow subworkflow_input {
   input {Int i Int j }
  call sub.subwf {input: i=i, j=j}
  meta {
    allowNestedInputs: true
  }
}
