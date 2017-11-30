workflow flatten {

  Array[Array[Int]] identity = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
  ]

  Array[Map[Int, String]] list_of_maps = [
    { 1: "one", 2: "two" },
    { 11: "eleven", 22: "twenty-two" }
  ]

  output {
    Array[Int] flattened_identity = flatten(identity)
    Array[Pair[Int, String]] flattened_map = flatten(list_of_maps)
  }
}
