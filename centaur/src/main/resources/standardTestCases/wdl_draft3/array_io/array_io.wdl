version 1.0

workflow array_io {
  Array[Int] rs = range(3)

  scatter (r in rs) {
    call mk_file { input: index = r }
  }

  Array[String] strings = ["str1", "str2", "str3"]
  Array[File] files = mk_file.out

  call serialize { input: strs = strings }
  call concat_files as concat { input: files = files }
  call count_lines { input: files = [concat.concatenated] }
  call count_lines as count_lines_array { input: files = files }

  output {
    Int array_count = count_lines_array.count
    Int single_count = count_lines.count
    String serialized = serialize.contents
    String concatenated = read_string(concat.concatenated)
  }
}

task concat_files {
  input {
    String? flags
    Array[File]+ files
  }
  command <<<
    cat ~{default = "-s" flags} ~{sep = " " files}
  >>>
  output {
    File concatenated = stdout()
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task count_lines {
  input {
    Array[File]+ files
  }

  command <<<
    cat ~{sep=" " files} | wc -l
  >>>
  output {
    Int count = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task serialize {
  input {
    Array[String] strs
  }
  command <<<
    cat ~{write_lines(strs)}
  >>>
  output {
    String contents = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task mk_file {
  input {
    Int index
  }
  command <<<
    echo "file-~{index}" > i
  >>>
  output {
    File out = "i"
  }
  runtime { docker: "ubuntu:latest" }
}
