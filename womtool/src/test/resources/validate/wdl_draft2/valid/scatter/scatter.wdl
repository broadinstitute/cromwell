workflow scatter_wf {
  scatter (x in range(2)) {
    Int i = 5
    call scatter_task { input: i = i }
  }

  call scatter_task as st2 { input: i = scatter_task.i_again[0] }

  output {
    Int nope = 5
    Array[Int] i_agains = scatter_task.i_again
    Int i_yet_again = st2.i_again
  }
}

task scatter_task {
  Int i
  command { echo ${i} }
  output {
    Int i_again = read_int(stdout())
  }
}
