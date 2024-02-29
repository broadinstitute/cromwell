version 1.0

import "scatter_chain.wdl" as chain

workflow w {
  call chain.scatter_chain
  output {
    Array[Int] incremented = scatter_chain.incremented
  }
}
