version development

workflow use_directory {
  call make_directory
  call read_from_directory { input: d = make_directory.d }

  output {
    Array[String] out = read_from_directory.contents
  }
}

task make_directory {
  command {
    mkdir foo
    echo "foo text" > foo/text
    echo "foo text2" > foo/text2
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Directory d = "foo/"
  }
}

task read_from_directory {
  input {
    Directory d
  }
  command {
    cat ~{d}/text
    cat ~{d}/text2
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Array[String] contents = read_lines(stdout())
  }
}
