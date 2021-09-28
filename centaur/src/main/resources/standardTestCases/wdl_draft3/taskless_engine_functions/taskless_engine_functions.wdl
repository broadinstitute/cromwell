version 1.0

workflow taskless_engine_functions {

  Array[Int] ints = [ 1, 2 ]

  Array[String] strings = ["a", "b"]

  String filepath = "gs://not/a/real/file.txt"

  Array[Array[Int]] matrix = [
    [1, 0],
    [1, 0]
  ]

  Array[Map[Int, String]] list_of_maps = [
    { 1: "one", 2: "two" },
    { 11: "eleven", 22: "twenty-two" }
  ]

  Float f = 1.024

  output {
    Array[Pair[Int, String]] int_cross_string = cross(ints, strings)
    Array[Array[Int]] transposed_matrix = transpose(matrix)

    Array[Int] flattened_matrix = flatten(matrix)
    Int matrix_length = length(matrix)
    Int flattened_matrix_length = length(flattened_matrix)
    Array[Pair[Int, String]] flattened_map = flatten(list_of_maps)

    String file_basename = basename(filepath)
    String file_basename_extensionless = basename(filepath, ".txt")

    Int f_floor = floor(f)
    Int f_ceiling = ceil(f)
    Int f_round = round(f)
  }
}
