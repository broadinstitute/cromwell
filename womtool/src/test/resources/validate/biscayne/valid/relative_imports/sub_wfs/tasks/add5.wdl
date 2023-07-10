version development-1.1

import "../../structs/my_struct.wdl"

task add5 {
  input {
    MyStruct x
  }
  command <<<
    echo $((5 + ~{x.a}))
  >>>
  output {
    MyStruct five_added = object { a: read_int(stdout()) }
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
