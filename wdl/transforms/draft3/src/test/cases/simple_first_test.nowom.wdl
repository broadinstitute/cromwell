version draft-3

task in_n_out {
  input {
    Int total
  }
  command { echo ${total} }
  output {
    Int out = read_int(stdout()) + 1
  }
}

workflow order {
  input {
    Int n = 4
  }
  call in_n_out { input: total = n }
}
