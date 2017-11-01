task file_outputs_from_input {
  String outputName1
  String outputName2
  String outputName3

  command {
    echo "foo" > ${outputName1}
    echo "bar" > ${outputName2}
    echo "baz" > ${outputName3}.txt
  }
  output {
    File foo = outputName1
    File bar = "${outputName2}"
    File baz = "${outputName3}.txt"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow file_outputs_from_input_wf {
  call file_outputs_from_input { input: outputName1 = "___something___", outputName2 = "___anything___", outputName3 = "___nothing___" }
  output {
    String foo = read_string(file_outputs_from_input.foo)
    String bar = read_string(file_outputs_from_input.bar)
    String baz = read_string(file_outputs_from_input.baz)
  }
}
