task cat {
  File asdf
  command {
    cat ${asdf}
  }
  output {
    File stuff = stdout()
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task grep {
  File in_file
  String pattern
  command {
    grep '${pattern}' ${in_file}
  }
  output {
    File blah = stdout()
    String x = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow composedenginefunctions {
  call cat
  call grep {
    input: in_file=cat.stuff
  }
  output {
    grep.x
  }
}
