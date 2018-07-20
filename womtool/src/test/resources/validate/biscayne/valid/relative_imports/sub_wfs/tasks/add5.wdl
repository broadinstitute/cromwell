version 1.0

import "../../structs/my_struct.wdl"

task add5 {
  input {
    MyStruct x
  }
  command <<<
    echo $((5 + ~{x}))
  >>>
  output {
    MyStruct five_added = object { a: read_int(stdout()) }
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
