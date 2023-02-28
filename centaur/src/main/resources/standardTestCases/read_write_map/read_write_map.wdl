task write_map {
  Map[String, String] file_to_name
  command {
    cat ${write_map(file_to_name)}
  }
  output {
    String contents = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task read_map {
  command <<<
    python <<CODE
    map = {'x': 500, 'y': 600, 'z': 700}
    print("\\n".join(["{}\\t{}".format(k,v) for k,v in map.items()]))
    CODE
  >>>
  output {
    Map[String, Int] out_map = read_map(stdout())
  }
  runtime {
    docker: "python:3.5.0"
  }
}

task assert_write_map_terminal_newline {
  Map[String, String] file_to_name
  command <<<
    # https://stackoverflow.com/a/25749716/21269164
    file_ends_with_newline() {
      [[ $(tail -c1 "$1" | wc -l) -gt 0 ]]
    }

    if ! file_ends_with_newline ${write_map(file_to_name)}
    then
      echo >&2 "Error: write_map() should write a file whose last character is a newline"
      exit 1
    fi
  >>>
  runtime {
    docker: "python:3.5.0"
  }
}

workflow wf {
  Map[String, String] map = {"f1": "alice", "f2": "bob", "f3": "chuck"}
  call write_map {input: file_to_name = map}
  call read_map
  call assert_write_map_terminal_newline {input: file_to_name = map}
  output {
     Map[String, Int] out = read_map.out_map
     String contents = write_map.contents
     }
}
