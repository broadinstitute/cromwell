version draft-3

workflow simple_scatter {
  Array[Int] indices = [1,2,3]

  scatter(i in indices) {
    Int j = i + 10
  }

  output {
    Array[Int] js = j
  }
}
