version 1.0

workflow nested_scatter {

  Array[Int] indices = [1,2,3]
  Int y = 55

  scatter(a in indices) {
    scatter(b in indices) {
      Int x = a + b
      scatter(c in indices) {
        Int j = a + b + c + x
      }
      scatter(d in j) {
        Int k = d + y
      }
    }
  }

  scatter (i in indices) {
    String string_in_scatter = "hello " + i
    if (i > 1) {
      String string_in_if_in_scatter = string_in_scatter
    }
  }

  output {
    Array[Array[Array[Int]]] ks = k
    Array[String?] strings_in_if_in_scatter = string_in_if_in_scatter
  }
}
