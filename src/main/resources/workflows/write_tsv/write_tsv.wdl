task tsv {
  Array[Array[String]] s = [["a", "b"],["c", "d"]]
  command {
    cat ${write_tsv(s)}
  }
  runtime {docker: "ubuntu:latest"}
  output {
    String o = read_string(stdout())
  }
}

workflow w {
  call tsv
}
