version 1.0

workflow two_level_scatter {

  Array[Int] indices = [1,2,3]

  scatter(a in indices) {
    scatter(b in indices) {
      Int x = a + b
    }
  }
}
