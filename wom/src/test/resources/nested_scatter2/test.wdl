task inner {
  Int i
  command { echo ${i} }
  output {
    Int out = read_int(stdout()) + 1
  }
}

task outer {
  Array[Array[Int]] matrix
  command { cat ${write_tsv(matrix)} }
  output {
    String tsv = read_string(stdout())
  }
}

workflow w {
  Array[Array[Int]] array = [[0,1,2],[3,4,5],[6,7,8]]

  scatter(i in array) {
    scatter(j in i) {
      call inner {input: i=j}
    }
  }

  # inner.out should be [[1,2,3],[4,5,6],[7,8,9]] (???)
  call outer {input: matrix=inner.out}
}
