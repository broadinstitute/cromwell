version 1.0

workflow biscayne_functions {
  Map[String, Int] abc = { "a": 1, "b": 2, "c": 3}

  output {
    # Oops! This is *not* a WDL 1.0 function!
    Array[String] abc_keys = keys(abc)
  }
}
