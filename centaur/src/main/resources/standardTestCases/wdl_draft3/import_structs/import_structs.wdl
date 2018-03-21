version draft-3

import "structs.wdl"

workflow import_structs {
  A a = object { i: 5, f: 5.5 }
  B b = object { a: a }

  output {
    B b_out = b
  }
}
