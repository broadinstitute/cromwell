version 1.2

struct Name {
  String first
  String last
}

workflow new_engine_functions {
  Array[Int] numbers = [1, 2, 3, 4, 5]
  Map[String, Int] scores = {"a": 1, "b": 2}
  Name name = Name { first: "Jane", last: "Doe" }
  Int? nothing = None

  output {
    # Operators
    Int int_power = 2 ** 10
    Float float_power = 2.0 ** 3
    Boolean none_equals_none = nothing == None

    # String functions
    String? found = find("abc123def456", "[0-9]+")
    String? not_found = find("abcdef", "[0-9]+")
    Boolean does_match = matches("abc123", "^a.+3$")
    String swapped = sub("hello world", "([a-z]+) ([a-z]+)", "\\2 \\1")

    # File functions
    String joined = join_paths("/a/b", ["c", "d.txt"])

    # Array functions
    Boolean has_three = contains(numbers, 3)
    Array[Array[Int]] chunks = chunk(numbers, 2)
    Int first_or_default = select_first([nothing], 42)

    # Map, struct and object functions
    Boolean has_a = contains_key(scores, "a")
    Boolean has_c = contains_key(scores, "c")
    Array[Int] score_values = values(scores)
    Array[String] name_keys = keys(name)
    Int scores_length = length(scores)
    Int string_length = length("hello")
  }
}
