version development

workflow directory_type {
  input {
    String text2loc = "text2"
  }
  call make_directory { input: text2loc = text2loc }
  call read_from_directory { input: d = make_directory.d, text2loc = text2loc }

  output {
    Array[String] out = read_from_directory.contents
  }
}

task make_directory {
  input {
    String text2loc
  }
  String text2dir = sub("foo/~{text2loc}", "/[^/]*$", "")
  command {
    mkdir foo
    mkdir -p ~{text2dir}
    echo "foo text" > foo/text
    echo "foo text2" > foo/~{text2loc}
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
    String text2loc
    Directory d
  }
  command {
    cat ~{d}/text
    cat ~{d}/~{text2loc}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Array[String] contents = read_lines(stdout())
  }
}
