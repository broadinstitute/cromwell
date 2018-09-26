version 1.0

struct MiniFoo {
  Int one
}

struct FooStruct {
  MiniFoo field1
  Map[String, String] field2
}

workflow draft3_read_and_coerce_json {
  FooStruct my_foo = object {
    field1: object { one: 1 },
    field2: {
      "field1": "value1",
      "field2": "value2"
    }
  }

  call read_write_all { input: my_foo = my_foo }

  output {
    String from_coerced_map = read_write_all.coerced_to_map.field2["field1"]
    String from_coerced_foo = read_write_all.coerced_to_foo.field2["field2"]
  }
}

task read_write_all {
  input {
    FooStruct my_foo
  }

  File written_json = write_json(my_foo)

  command <<<
    cp ~{written_json} copied_written.json
    cp ~{written_json} copied_written_2.json
  >>>

  output {
    Map[String, Map[String, String]] coerced_to_map = read_json("copied_written.json")
    FooStruct coerced_to_foo = read_json("copied_written_2.json")
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
