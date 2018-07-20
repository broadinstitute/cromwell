version 1.0

import "sub_wfs/foo.wdl"

workflow relative_imports {
  call foo.foo_wf

  output {
    Int result = foo_wf.unpacked
  }
}
