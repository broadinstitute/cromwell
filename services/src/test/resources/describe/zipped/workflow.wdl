version 1.0

import "simple/simple.wdl" as simple_pkg

workflow use_simple {
  input {
    Int i
  }
  call simple_pkg.simple { input: i = i }
  output {
    Int result = simple_pkg.simple.result
  }
}