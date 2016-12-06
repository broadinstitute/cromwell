task prepare {
 command <<<
    python3 -c "print('one\ntwo\nthree\nfour')"
    sleep 2
  >>>
  output {
    Array[String] array = read_lines(stdout())
  }
  runtime {
    docker: "python:3.5.0"
  }
}

task analysis {
  String str
  command <<<
    python3 -c "print('_${str}_')" > a.txt
    sleep 2
  >>>
  output {
    File out = "a.txt"
  }
  runtime {
    docker: "python:3.5.0"
  }
}

task gather {
  Array[File] array
  command <<<
    cat ${sep=' ' array}
    sleep 2
  >>>
  output {
    String str = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow scattergather {
  call prepare
  scatter (x in prepare.array) {
    call analysis {input: str=x}
  }
  call gather {input: array=analysis.out}
  output {
    gather.str
  }
}
