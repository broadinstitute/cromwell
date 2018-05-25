version 1.0

workflow subwf {
  Int i = 5

  input {
    Int j = i + 1000
  }

  output {
    Int k = j
  }
}
