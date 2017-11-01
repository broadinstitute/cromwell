task golden_pie {
  Float pi = 3.1415926
  Float tau = pi + pi

  command {
    echo 1.6180339887
    echo ${tau} 1>&2
  }

  output {
    Float Au = read_float(stdout())
    Float doubleAu = Au + Au
    Float tauValue = read_float(stderr())
  }

  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf {
  call golden_pie
}
