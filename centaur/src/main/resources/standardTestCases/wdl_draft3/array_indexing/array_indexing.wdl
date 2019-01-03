version 1.0

struct Indexable {
  Int int_field
}

workflow array_indexing {
  Array[Int] to_index = range(35)

  Indexable indexable_struct = object { int_field: 55 }

  Map[String, Int] string_indexed_map = { "int_value": 101 }

  Map[Int, String] int_indexed_map = { 101: "int_value" }

  output {
    # to_index is an identity function so this should dereference to 25: (5 + (50-30)):
    Int array_deref = to_index[25]
    Int struct_deref = indexable_struct["int_field"]
    Int string_map_deref = string_indexed_map["int_value"]
    String int_map_deref = int_indexed_map[101]
  }
}
