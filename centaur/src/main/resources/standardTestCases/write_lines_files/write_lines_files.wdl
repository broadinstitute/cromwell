task make_file {
  Int index

  command {
    echo contents_${index} > out
  }

  output {
    File out = "out"
  }
  runtime { docker:"ubuntu:latest" }
}

task switcho_reverso {
  Array[File] files

  command {
    F="${write_lines(files)}"
    for f in $(tac $F)
    do
      cat $f
    done
  }

  output {
    Array[String] out = read_lines(stdout())
  }
  runtime {
    docker:"ubuntu:latest"
    failOnStderr: true
  }
}

workflow write_lines_files {
  Array[Int] is = range(5)
  scatter(i in is) {
    call make_file { input: index = i }
  }

  call switcho_reverso { input: files = make_file.out }

  output {
    Array[String] reversed_contents = switcho_reverso.out
  }
}
