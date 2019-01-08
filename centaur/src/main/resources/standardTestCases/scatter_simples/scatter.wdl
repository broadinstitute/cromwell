version 1.0

import "bag_of_tasks.wdl" as bag

workflow w {
  call bag.A
  scatter (item in A.A_out) {
    call bag.B {input: B_in = item}
    call bag.C {input: C_in = B.B_out}
    call bag.E
  }
  call bag.D {input: D_in = B.B_out}
}
