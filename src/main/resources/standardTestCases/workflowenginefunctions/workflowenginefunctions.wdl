task printArray {
  Array[Int] array

  command {
    echo "${sep='"; echo "' array}"
    sleep 2
  }
  output {
    String o = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow workflowEngineFunctions {
  Array[Int] some_ints = range(10)
  call printArray {input: array=some_ints}
}
