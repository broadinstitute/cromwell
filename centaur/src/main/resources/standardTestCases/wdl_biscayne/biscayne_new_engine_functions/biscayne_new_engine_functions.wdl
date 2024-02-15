version development-1.1

workflow biscayne_new_engine_functions {

  meta {
    description: "This test makes sure that these functions work in a real workflow"
    functions_under_test: [ "keys", "as_map", "as_pairs", "collect_by_key", "suffix" ]
  }

  Map[String, Int] x_map_in = {"a": 1, "b": 2, "c": 3}
  Map[String,Pair[File,File]] y_map_in = { "a": ("a.bam", "a.bai"), "b": ("b.bam", "b.bai") }

  Array[Pair[String,Int]] x_pairs_in = [("a", 1), ("b", 2), ("c", 3)]
  Array[Pair[String,Pair[File,File]]] y_pairs_in = [("a", ("a.bam", "a.bai")), ("b", ("b.bam", "b.bai"))]

  Array[Pair[String,Int]] z_pairs_in = [("a", 1), ("b", 2), ("a", 3)]

  Array[String] some_strings = ["aaa", "bbb", "ccc"]

  Int smallestInt = 1
  Float smallFloat = 2.718
  Float bigFloat = 3.141
  Int biggestInt = 10

  Int maxInt = 2147483647
  # max float... near enough:
  Float maxFloat = 179769313000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.0

  output {

    # keys(), as_map(), as_pairs(), collect_by_key():
    # ==============================================
    Array[String] x_map_keys_out = keys(x_map_in) # ["a", "b", "c"]
    Array[String] y_map_keys_out = keys(y_map_in) # ["a", "b"]

    Array[Pair[String,Int]] x_as_pairs_out = as_pairs(x_map_in) # [("a", 1), ("b", 2), ("c", 3)]
    Array[Pair[String,Pair[File,File]]] y_as_pairs_out = as_pairs(y_map_in) # [("a", ("a.bam", "a.bai")), ("b", ("b.bam", "b.bai"))]

    Map[String,Int] x_as_map_out = as_map(x_pairs_in) # {"a": 1, "b": 2, "c": 3}
    Map[String,Pair[File,File]] y_as_map_out = as_map(y_pairs_in) # {"a": ("a.bam, "a.bai"), "b": ("b.bam", "b.bai")}

    Map[String,Array[Int]] x_collect_by_key_out = collect_by_key(x_pairs_in) # {"a": [1], "b": [2], "c": [3]}
    Map[String,Array[Int]] z_collect_by_key_out = collect_by_key(z_pairs_in) # {"a": [1, 3], "b": [2]}

    # min(), max():
    # =================================================
    Int smallerInt = min(biggestInt, smallestInt) # 1
    Int biggerInt = max(smallestInt, biggestInt) # 10
    Float smallIntFloatComparison = min(smallestInt, smallFloat) # 1.0
    Float bigIntFloatComparison = max(bigFloat, biggestInt) # 10.0
    Float minMaxIntFloatComposition = min(max(biggestInt, smallFloat), smallestInt) # 1.0
    Float maxIntVsMaxFloat = max(maxInt, maxFloat)

    # suffix():
    # =================================================
    Array[String] with_suffixes = suffix("S", some_strings)
  }
}

