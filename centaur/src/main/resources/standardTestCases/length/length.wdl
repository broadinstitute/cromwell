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

workflow length {
  Array[String] empty = []

  call lens

  output {
    Int someStringsLength = length(lens.someStrings)
    Int someIntsLength = length(lens.someInts)
    Int emptyLength = length(empty)
  }
}
