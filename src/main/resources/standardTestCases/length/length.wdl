task lens {
  command {
    echo whatevs
  }
  output {
    Array[String] someStrings = ["foo", "bar", "baz"]
    Array[Int] someInts = [1, 2, 3]
    Array[File] noFiles = glob("*.no.such.file.txt")
  }
}

workflow length {
  call lens
  output {
    Int someStringsLength = length(lens.someStrings)
    Int someIntsLength = length(lens.someInts)
    Int noFilesLength = length(lens.noFiles)
  }
}
