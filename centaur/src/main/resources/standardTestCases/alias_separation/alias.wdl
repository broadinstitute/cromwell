version 1.0

workflow alias_separation {
  call t as t1 { input: str="t1" }
  call t as t2 { input: str="t2" }

  output {
    Array[String] t1out = t1.out
    Array[String] t2out = t2.out
  }
}

task t {
  input {
    String str
  }

  command {
    echo ${str} >> outfile
  }

  output {
    Array[String] out = read_lines("outfile")
  }

  runtime {
    docker: "amazonlinux:latest"
  }
}
