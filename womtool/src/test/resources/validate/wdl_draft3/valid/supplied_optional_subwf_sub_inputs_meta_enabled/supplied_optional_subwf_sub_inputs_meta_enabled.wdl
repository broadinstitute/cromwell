version 1.0

import "import_me.wdl"

workflow supplied_optional_subwf_sub_inputs {
  input {
    Int y
  }
  # Shouldn't (strictly) be able to call this sub-workflow because the inputs are not passed through
  call import_me.sub_wf {input: y=y}
  meta {
    allowNestedInputs: true
  }
}
