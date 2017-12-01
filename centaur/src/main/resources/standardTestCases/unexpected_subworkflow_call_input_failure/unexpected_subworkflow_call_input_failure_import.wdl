workflow subwf {
  Int i = 5
  Int j = i + 1000

  output {
    Int k = j
  }
}
