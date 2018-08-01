workflow bad_autobox {
  Int i = 5

  # Bad autobox!
  Array[Int] bad_is = i

  # Good manual-box!
  Array[Int] good_is = [ i ]
}
