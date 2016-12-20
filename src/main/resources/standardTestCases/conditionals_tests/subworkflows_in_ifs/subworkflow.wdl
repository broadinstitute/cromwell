task increment {
  Int i
  command {
    echo $(( ${i} + 1 ))

    sleep 2
  }
  output {
    Int j = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow subwf {
  Array[Int] is
  scatter (i in is) {
    call increment { input: i = i }
  }
  output {
    Array[Int] js = increment.j
  }
}