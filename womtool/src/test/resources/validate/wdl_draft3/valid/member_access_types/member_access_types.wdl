version 1.0

workflow member_access {
  # Object access:
  Object myObj = object { an_int: 5 }
  if (myObj.an_int == 10) {
    Boolean a = true
  }

  # Array access:
  Array[Int] is = [0, 1]
  if (is[1] == 1) {
    Boolean b = true
  }

  # Map access:
  Map[String, Int] msb = { "an_int": 55 }
  if (msb["an_int"] == 55) {
    Boolean c = true
  }

  # Array[Pair] access:
  Array[Pair[Boolean, Boolean]] apbb = [ (true, false), (false, true) ]
  if (apbb[0].left && apbb[1].right) {
    Boolean d = true
  }

  Object wrapped_array = object { an_array: [0, 1, 2] }
  scatter (x in wrapped_array.an_array) {
    Boolean e = true
  }

  output {
    Boolean a_out = select_first([a, false])
    Boolean b_out = select_first([b, false])
    Boolean c_out = select_first([c, false])
    Boolean d_out = select_first([d, false])
    Boolean e_out = e[0]
  }
}
