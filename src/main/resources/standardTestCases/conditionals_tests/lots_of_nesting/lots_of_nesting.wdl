import "ifs_in_scatters.wdl" as iis

workflow lots_of_nesting {

  if (true) {
    call iis.ifs_in_scatters as ifs_in_scatters_in_truth

    Array[Int?] allMirrors = ifs_in_scatters_in_truth.mirrors # [ null, 2, null, 4, null ]
    Array[Int] goodMirrors = select_all(allMirrors) # [2, 4]

    Array[Pair[Int, Int]] mxm = cross(goodMirrors, goodMirrors) # [ (2, 2) (2, 4), (4, 2), (4, 4) ]

    scatter (m in mxm) {
      Int lhs = m.left
      if (lhs * m.right == 8) {
        Int rhs = m.right
        Pair[Int, Int] hooroo = ( lhs, rhs )
      }
      Pair[Int, Int?] hooray = ( lhs, rhs )
    }
  }

  if (false) {
    call iis.ifs_in_scatters as ifs_in_scatters_in_falsehood
  }


  output {
    Array[Int?]? a = allMirrors
    Array[Int]? b = goodMirrors
    Array[Pair[Int, Int]]? c = mxm
    Array[Pair[Int, Int?]]? d = hooray
    Array[Pair[Int, Int]?]? e = hooroo
  }
}
