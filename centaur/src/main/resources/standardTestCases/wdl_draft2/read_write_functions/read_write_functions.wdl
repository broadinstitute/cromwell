workflow d2_read_write_functions {
  Object my_foo = object {
    field1: "1",
    field2: 2
  }

  call read_write_all { input: my_foo = my_foo }

  output {
    Array[String] lines = read_write_all.lines
    Array[Array[String]] tsv = read_write_all.tsv
    Object json = read_write_all.json
  }
}

task read_write_all {
  Object my_foo

  String line3 = "line3"

  Array[String] lines_to_write = ["line1", "line2", "line3"]
  Array[Object] tsv_to_write = [ object { line1: "line one", line2: "line two", line3: "line three" } ]

  File written_lines = write_lines(lines_to_write)
  File written_tsv = write_tsv(tsv_to_write)
  File written_json = write_json(my_foo)

  command <<<
    mv ${written_lines} written_lines2
    mv ${written_tsv} written_tsv2
    mv ${written_json} written_json2
  >>>

  output {
    Array[String] lines = read_lines("written_lines2")
    Array[Array[String]] tsv = read_tsv("written_tsv2")
    Object json = read_json("written_json2")
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
