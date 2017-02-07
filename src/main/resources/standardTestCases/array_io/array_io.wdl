task concat_files {
  String? flags
  Array[File]+ files
  command {
    cat ${default = "-s" flags} ${sep = " " files} && sleep 2
  }
  output {
    File concatenated = stdout()
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task count_lines {
  Array[File]+ files
  command {
    cat ${sep=" " files} | wc -l && sleep 2
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task serialize {
  Array[String] strs
  command {
    cat ${write_lines(strs)} && sleep 2
  }
  output {
    String contents = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task mk_file {
  Int index
  command { echo "file-${index}" > i && sleep 2 }
  output { File out = "i" }
  runtime { docker: "ubuntu:latest" }
}

workflow array_io {
  Array[Int] rs = range(3)
  scatter (r in rs) {
    call mk_file { input: index = r }
  }
  Array[String] strings = ["str1", "str2", "str3"]
  call serialize {
    input: strs = strings
  }
  call concat_files as concat {
    input: files = mk_file.out
  }
  call count_lines {
    input: files = [concat.concatenated]
  }
  call count_lines as count_lines_array {
    input: files = mk_file.out
  }

  output {
    Int array_count = count_lines_array.count
    Int single_count = count_lines.count
    String serialized = serialize.contents
    String concatenated = read_string(concat.concatenated)
  }
}
