version 1.0

workflow lots_of_nesting {

  Boolean b0 = true
  Boolean b1 = true
  Boolean b2 = true

  Array[Int] i0s = range(2)
  Array[Int] i1s = range(2)
  Array[Int] i2s = range(2)

  String s0 = "hello"
  String s1 = "world"

  if (b0) {
    scatter(i0 in i0s) {
      if (b1) {
        scatter(i1 in i1s) {
          if (b2) {
            scatter(i2 in i2s) {
              String s = s0 + s1
            }
          }
        }
      }
    }
  }

  output {
    Array[Array[Array[String]?]?]? s_out = s
  }
}
