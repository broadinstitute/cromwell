workflow read_write_json {

  Int? supplied_optional_int = 1
  Int? not_supplied_optional_int

  Object json_object = object {
    str: "hi",
    emptystr: "",
    int: 57,
    float: 27.5,
    pair: (5, "hello"),
    array: ["a", "b", "c"],
    obj: object {
      inner_str: "inner hello",
      inner_float: 28.5
    },
    opt_int_1: supplied_optional_int,
    opt_int_2: not_supplied_optional_int
  }

  call make_some_json
  call round_trip { input: to_jsonify = json_object }

  Object round_tripped_inner_obj = round_trip.round_tripped.obj

  if (round_trip.round_tripped.emptystr == "" &&
    round_tripped_inner_obj.inner_float == make_some_json.output_json.float + round_trip.round_tripped.opt_int_1) {
      call success { input: actual = round_tripped_inner_obj.inner_float, expected = make_some_json.output_json.float }
  }

  output {
    String? result = success.result
  }
}

task round_trip {
  Object to_jsonify
  File jsonified = write_json(to_jsonify)
  command {
    cp ${jsonified} output.json
  }
  output {
    Object round_tripped = read_json("output.json")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task make_some_json {
  command <<<
    echo '{'
    echo '  "str": "hi",'
    echo '  "emptystr": "",'
    echo '  "int": 57,'
    echo '  "float": 27.5,'
    echo '  "pair": {'
    echo '    "left": 5,'
    echo '    "right": "hello"'
    echo '  },'
    echo '  "array": ["a", "b", "c"]'
    echo '}'
  >>>
  output {
    Object output_json = read_json(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task success {
  Float actual
  Float expected
  command {
    echo "${actual} was equal to ${expected} plus 1"
  }
  output {
    String result = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
