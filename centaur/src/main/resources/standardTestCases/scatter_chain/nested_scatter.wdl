version 1.0

import "scatter_chain.wdl" as chain

workflow w {
  call chain
  output {
    Array[Int] incremented = chain.incremented
  }
}
