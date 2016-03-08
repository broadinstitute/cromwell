task printArray {
  Array[String] array

  command {
    echo "${sep='"; echo "' array}"
  }
  output {
    String o = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow workspaceEngineFunctions {
  File filein
  Array[String] in_file = read_lines(filein)
  call printArray {input: array=in_file}
}
