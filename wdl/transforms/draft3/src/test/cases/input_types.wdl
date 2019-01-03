version 1.0

workflow input_types {
  input {
    # All the primitive types:
    Int i
    String s
    Float float
    Boolean b
    File file
    Object o

    # Simple compound types:
    Int? maybe_i
    Array[String] array_s
    Map[Int, String] map_is

    # Check complex nesting:
    Array[Pair[String, Int]?] lotsa_nesting_array
  }
}
