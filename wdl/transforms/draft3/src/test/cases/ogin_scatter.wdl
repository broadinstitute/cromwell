version 1.0

workflow ogin_scatter {

  Array[Int] indices = [1,2,3]
  Int ogin_me = 10

  scatter(i in indices) {
    Int j = i + ogin_me
  }

  output {
    Array[Int] js = j
  }
}
