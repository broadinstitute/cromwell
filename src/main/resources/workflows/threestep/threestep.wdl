task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task cgrep {
  String pattern
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow three_step {
  call ps
  call cgrep { input: in_file=ps.procs }
  call wc { input: in_file=ps.procs }
  output {
    wc.count
  }
}
