task grep {
  command {
    grep ${pattern} ${flags?} ${file file_name}
  }
  output {
    file out = "stdout"
  }
  runtime {
    memory: "2MB"
    cores: 1
    disk: "3MB"
  }
}

task wc {
  command {
    wc -l ${sep=' ' file files+} | tail -1 | cut -d' ' -f 2
  }
  output {
    int count = read_int("stdout")
  }
}

workflow scatter_gather_grep_wc {
  array[file] input_files
  scatter(f in input_files) {
    call grep {
      input: file_name = f
    }
  }
  call wc {
    input: files = grep.out
  }
}


