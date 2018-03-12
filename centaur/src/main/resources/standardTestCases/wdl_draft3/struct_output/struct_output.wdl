version draft-3

struct FooStruct {
  Int simple
  Pair[Array[Int], Map[String, Boolean]] complex
}

workflow struct_output {
  output {
    FooStruct myFoo = object {
      simple: 5,
      complex: ([5], { "t": true })
    }
  }
}
