version 1.0

workflow scatter_var_member_access {
  Array[Pair[Int, Int]] pairs = [(1, 2), (3, 4), (5, 6)]
  scatter (p in pairs) {
    Int x = p.left
  }
}
