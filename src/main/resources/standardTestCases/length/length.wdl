task lens {
  command {
    echo whatevs
  }
  output {
    Array[String] someStrings = ["foo", "bar", "baz"]
    Array[Int] someInts = [1, 2, 3]
  }
  runtime { docker: "ubuntu:latest" }
}

task void {
  command {}
  output {
    String out = ""
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow length {
  # Hack to create an empty Array
  if (false) { call void }
  Array[String] empty = select_all([void.out])

  call lens

  output {
    Int someStringsLength = length(lens.someStrings)
    Int someIntsLength = length(lens.someInts)
    Int emptyLength = length(empty)
  }
}
