task a2f {
  Array[Array[String]] strings=[["a","b"],["c","d"]]
  File tsv = write_tsv(strings)

  command {
    cat ${write_tsv(strings)} > f0
    cat ${tsv} > f1
    sleep 2
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    Array[Array[String]] out0 = read_tsv("f0")
    Array[Array[String]] out1 = read_tsv("f1")
    String out2 = read_string("f0")
  }
}

workflow write_tsv {
  call a2f
}
