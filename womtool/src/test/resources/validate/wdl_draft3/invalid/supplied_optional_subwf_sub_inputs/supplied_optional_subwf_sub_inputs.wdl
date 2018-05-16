version 1.0

import "import_me.wdl"

workflow supplied_optional_subwf_sub_inputs {
  # Shouldn't (strictly) be able to call this sub-workflow because the inputs are not passed through
  call import_me.sub_wf
}
