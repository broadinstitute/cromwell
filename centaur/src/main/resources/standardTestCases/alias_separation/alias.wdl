version 1.0

workflow alias_separation {

  scatter (x in range(3)) {
    call t as t1 { input: str="t1" }
    call t as t2 { input: str="t2" }
  }

  output {
    Array[Array[String]] t1out = t1.out
    Array[Array[String]] t2out = t2.out
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
    File outfile = "outfile"
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
