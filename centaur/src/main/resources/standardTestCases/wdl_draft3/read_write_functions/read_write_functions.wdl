version 1.0

struct FooStruct {
  String field1
  Int field2
  String field3
}

workflow read_write_functions {
  FooStruct my_foo = object {
    field1: "1",
    field2: 2,
    field3: ""
  }

  call read_write_all { input: my_foo = my_foo }

  output {
    Array[String] lines = read_write_all.lines
    Array[Array[String]] tsv = read_write_all.tsv
    Map[String, String] map = read_write_all.map
    Object my_object = read_write_all.my_object
    Array[Object] objects = read_write_all.objects
    FooStruct json = read_write_all.json
  }
}

task read_write_all {
  input {
    FooStruct my_foo
  }

  String line3 = "line3"

  File written_lines = write_lines(["line1", "line2", "line3"])
  File written_tsv = write_tsv( [[ "line1", "line one"], ["line2", "line two"], [ "line3", "line three" ]] )
  File written_map = write_map( {"key1": "value1", "key2": "value2", "key3": "value3"} )
  File written_object = write_object(object { line1: "line one", line2: "line two", line3: "line three" })
  File written_objects = write_objects( [object { line1: "line one", line2: "line two", line3: "line three" },
                                         object { line1: "line one (second obj)", line2: "line two (second obj)", line3: "line three (second obj)" }] )
  File written_json = write_json(my_foo)

  command <<<
    mv ~{written_lines} written_lines2
    mv ~{written_tsv} written_tsv2
    mv ~{written_map} written_map2
    mv ~{written_object} written_object2
    mv ~{written_objects} written_objects2
    mv ~{written_json} written_json2
  >>>

  output {
    Array[String] lines = read_lines("written_lines2")
    Array[Array[String]] tsv = read_tsv("written_tsv2")
    Map[String, String] map = read_map("written_map2")
    Object my_object = read_object("written_object2")
    Array[Object] objects = read_objects("written_objects2")
    FooStruct json = read_json("written_json2")
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
