version draft-3

import "structs.wdl"

workflow import_structs {
  A a = object { i: 5, f: 5.5 }
  B b = object { a: a }
  StructCollision sc = object { i: 5 }
  StructCollision2 sc2 = object { s: "hi" }

  output {
    B b_out = b
    String sc_out = sc.i
    String sc2_out = sc2.s
  }
}
