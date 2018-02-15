version draft-3

workflow input_types {
  input {
    # All the primitive types:
    Int i
    String s
    Float f
    Boolean b
    File f
    Object o

    # Simple compound types:
    Int? maybe_i
    Array[String] array_s
    Map[Int, String] map_is

    # Check complex nesting:
    Array[Pair[String, Int]?] lotsa_nesting_array
  }
}
