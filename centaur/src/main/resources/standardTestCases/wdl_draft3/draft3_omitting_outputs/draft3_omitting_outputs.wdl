version 1.0

import "import_me.wdl" as blah

workflow draft3_omitting_outputs {
  call blah.foo

  if (true) {
    scatter (x in range(1)) {
      call blah.foo_in_a_box as foo2
    }
  }
}
