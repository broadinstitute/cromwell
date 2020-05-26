version 1.0

import "import_me.wdl"

workflow supplied_optional_subwf_sub_inputs {
  input {Int y}
  call import_me.sub_wf {input: y=y}
}
