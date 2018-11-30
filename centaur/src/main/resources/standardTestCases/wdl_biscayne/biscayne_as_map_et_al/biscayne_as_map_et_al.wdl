version development

workflow biscayne_as_map_et_al {

  meta {
    description: "This test makes sure that these functions work in a real workflow"
    functions_under_test: [ "keys", "as_map", "as_pairs", "collect_by_key" ]
  }

  Map[String, Int] x_map_in = {"a": 1, "b": 2, "c": 3}
  Map[String,Pair[File,File]] y_map_in = { "a": ("a.bam", "a.bai"), "b": ("b.bam", "b.bai") }

  Array[Pair[String,Int]] x_pairs_in = [("a", 1), ("b", 2), ("c", 3)]
  Array[Pair[String,Pair[File,File]]] y_pairs_in = [("a", ("a.bam", "a.bai")), ("b", ("b.bam", "b.bai"))]

  Array[Pair[String,Int]] z_pairs_in = [("a", 1), ("b", 2), ("a", 3)]


  output {
    Array[String] x_map_keys_out = keys(x_map_in) # ["a", "b", "c"]
    Array[String] y_map_keys_out = keys(y_map_in) # ["a", "b"]

    Array[Pair[String,Int]] x_as_pairs_out = as_pairs(x_map_in) # [("a", 1), ("b", 2), ("c", 3)]
    Array[Pair[String,Pair[File,File]]] y_as_pairs_out = as_pairs(y_map_in) # [("a", ("a.bam", "a.bai")), ("b", ("b.bam", "b.bai"))]

    Map[String,Int] x_as_map_out = as_map(x_pairs_in) # {"a": 1, "b": 2, "c": 3}
    Map[String,Pair[File,File]] y_as_map_out = as_map(y_pairs_in) # {"a": ("a.bam, "a.bai"), "b": ("b.bam", "b.bai")}

    Map[String,Array[Int]] x_collect_by_key_out = collect_by_key(x_pairs_in) # {"a": [1], "b": [2], "c": [3]}
    Map[String,Array[Int]] z_collect_by_key_out = collect_by_key(z_pairs_in) # {"a": [1, 3], "b": [2]}
  }
}

