version development-1.1

import "../structs/my_struct.wdl"
import "tasks/add5.wdl" as a5

workflow foo_wf {
  call a5.add5 { input: x = object { a: 100 } }
  output {
    Int unpacked = add5.five_added.a
  }
}
