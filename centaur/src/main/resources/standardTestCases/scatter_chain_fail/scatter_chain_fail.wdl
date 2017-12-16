workflow scatter_chain_fail {
  Array[Int] is = [ 1, 2, 3, 4, 5 ]

  scatter (i in is) {
    call fail as i1 { input: i = i }
  }

  # This scatter will never be runnable because i1.0 will fail
  scatter (j in i1.o) {
    call fail as i2 { input: i = j }
  }

  output {
    Array[Int] workflow_o = i2.o
  }
}


task fail {
  Int i

  command {
    echo $(( ${i} + 1 ))
    exit 1
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Int o = read_int(stdout())
  }
}
