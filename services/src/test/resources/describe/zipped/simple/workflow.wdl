version 1.0

import "simple/simple.wdl" as simple_pkg

workflow use_simple {
  input {
    Int i = 1
  }
  call simple_pkg.simple { input: i = i }
  output {
    Int result = simple.result
  }
}