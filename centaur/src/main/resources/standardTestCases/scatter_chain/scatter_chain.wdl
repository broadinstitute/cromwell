version 1.0

task increment {
  input {
    Int i
  }

  command {
    echo $(( ${i} + 1 ))
    exit 1
  }
  runtime {
    docker: "shuubuntu:oldest"
  }
  output {
    Int o = read_int(stdout())
  }
}

task sum {
  input {
    Int i
    Int j
  }
  command {
    echo $(( ${i} + ${j} ))
    exit 1
  }
  runtime {
    docker: "shuubuntu:oldest"
  }
  output {
    Int o = read_int(stdout())
  }
}


workflow scatter_chain {
  Array[Int] is = [ 1, 2, 3, 4, 5 ]

  scatter (i in is) {
    call increment as i1 { input: i = i }
  }

  scatter (j in i1.o) {
    call increment as i2 { input: i = j }
  }

  Array[Pair[Int, Int]] zipped = zip(i1.o, i2.o)

  scatter (p in zipped) {
    call sum { input: i = p.left, j = p.right }
    call increment as i3 { input: i = p.left }
  }

  output {
    Array[Int] incremented = i3.o
    Array[Int] summed = sum.o
  }
}