version draft-3

task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}

task cgrep {
  input {
    String pattern
    File in_file
  }

  command {
    grep '~{pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

task wc {
  input {
    File in_file
  }
  command {
    cat ~{in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

workflow three_step {
  call ps
  call cgrep {
    input: in_file = ps.procs
  }
  call wc {
    input: in_file = ps.procs
  }
}
