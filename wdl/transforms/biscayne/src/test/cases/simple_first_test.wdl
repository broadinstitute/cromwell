version development-1.1

task in_n_out {
  input {
    Int total
    String amount
  }
  command { echo ${total} }
  output {
    Int out = read_int(stdout()) + 1
  }
}

workflow order {
  input {
    Int n = 4
    String more = "more"
  }
  call in_n_out { input: total = n, amount = more }

  output {
    Int out = in_n_out.out
  }
}
