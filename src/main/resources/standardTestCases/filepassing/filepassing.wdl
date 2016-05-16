task a {
  File in_file
  String out_name = "out"

  command {
    cat ${in_file} > ${out_name}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    File out = "out"
    File out_interpolation = "${out_name}"
    String contents = read_string("${out_name}")
  }
}

workflow filepassing {
  File f

  call a {input: in_file=f}
  call a as b {input: in_file=a.out}
  output {
      b.contents
  }
}
