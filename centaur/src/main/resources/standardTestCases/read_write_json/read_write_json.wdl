workflow read_write_json {
  Object json_object = object {
    str: "hi",
    int: 57,
    float: 27.5,
    pair: (5, "hello"),
    array: ["a", "b", "c"]
  }

  call make_some_json
  call round_trip { input: to_jsonify = json_object }

  if (round_trip.round_tripped.float == make_some_json.output_json.float) {
    call success { input: actual = round_trip.round_tripped.float, expected = make_some_json.output_json.float }
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
}

task make_some_json {
  command <<<
    echo '{'
    echo '  "str": "hi",'
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
}

task success {
  Float actual
  Float expected
  command {
    echo "${actual} was equal to ${expected}"
  }
  output {
    String result = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
