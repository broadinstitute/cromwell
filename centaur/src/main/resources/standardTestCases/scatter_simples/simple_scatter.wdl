version 1.0

task echo_int {
  input {
    Int int
  }
  command {
    echo ~{int}
  }
  output {
    Int out = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow scatter0 {
  Array[Int] ints = [1,2,3,4,5]
  call echo_int as outside_scatter {input: int = 8000}
  scatter(i in ints) {
    call echo_int as inside_scatter {
      input: int = i
    }
  }
}
