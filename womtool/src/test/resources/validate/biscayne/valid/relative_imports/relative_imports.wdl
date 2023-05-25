version development-1.1

import "sub_wfs/foo.wdl"

workflow relative_imports {
  call foo.foo_wf

  output {
    Int result = foo_wf.unpacked
  }
}
