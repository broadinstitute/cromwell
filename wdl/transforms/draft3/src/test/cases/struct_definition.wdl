version 1.0

struct FooStruct {
  Int simple
  Pair[Array[Int], Map[String, Boolean]] complex
}

workflow struct_definition {
  output {
    FooStruct myFoo = object {
      simple: 5,
      complex: ([5], { "t": true })
    }
  }
}
