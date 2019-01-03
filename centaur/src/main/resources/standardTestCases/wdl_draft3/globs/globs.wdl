version 1.0

workflow globs {
  call glob_output
  scatter(f in glob_output.generateds) {
    call int_reader { input: f = f }
  }
  call sum as summer { input: ints = int_reader.int }
  output {
    Int sum = summer.sum
  }
}

task sum {
  input {
    Array[Int]+ ints
  }
  command <<<
    echo $(( 0 + ~{sep=" + " ints} ))
  >>>
  output {
    Int sum = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task int_reader {
  input {
    File f
  }
  command {
    cat ~{f}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Int int = read_int(stdout())
  }
}

task glob_output {
  command <<<
    for i in `seq 0 5`
    do
      echo $i > text_${RANDOM}
    done
  >>>

  output {
    Array[File]+ generateds = glob("text_*")
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
