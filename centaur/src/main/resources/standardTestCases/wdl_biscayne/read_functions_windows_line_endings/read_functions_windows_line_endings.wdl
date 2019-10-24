version development

struct JsonObj {
    String field1
    Int field2
}

workflow biscayne_read_functions_windows_line_endings {

  call read_map
  call read_lines
  call read_tsv
  call read_json

  output {
    Map[String, Int] map = read_map.map
    Array[String] lines = read_lines.lines
    Array[Array[String]] tsv = read_tsv.tsv
    JsonObj json = read_json.json
  }
}

task read_map {
  command <<<
    python <<CODE
    map = {'x': 500, 'y': 600, 'z': 700}
    print("\\r\\n".join(["{}\\t{}".format(k,v) for k,v in map.items()]))
    CODE
  >>>

  output {
    Map[String, Int] map = read_map(stdout())
  }

  runtime {
    docker: "python:3.5.0"
  }
}

task read_lines {
  command <<<
    python <<CODE
    print("line1\\r\\nline2\\r\\nline3\\nline4")
    CODE
  >>>

  output {
    Array[String] lines = read_lines(stdout())
  }

  runtime {
    docker: "python:3.5.0"
  }
}

task read_tsv {
  command <<<
    python <<CODE
    print("line1\\tline one\\r\\nline2\\tline two\\r\\nline3\\tline three")
    CODE
  >>>

  output {
    Array[Array[String]] tsv = read_tsv(stdout())
  }

  runtime {
    docker: "python:3.5.0"
  }
}

task read_json {
  command <<<
    python <<CODE
    print("{\\r\\n\\"field1\\":\\"1\\",\\r\\n\\"field2\\":2\\r\\n}")
    CODE
  >>>

  output {
    JsonObj json = read_json(stdout())
  }

  runtime {
    docker: "python:3.5.0"
  }
}
