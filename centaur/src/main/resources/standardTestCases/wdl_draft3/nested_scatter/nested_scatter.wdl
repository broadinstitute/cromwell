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

  output {
    Array[Array[Array[Int]]] ks = k
  }
}
