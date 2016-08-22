task A {
  command {
    python -c "print(321);exit(123)"
  }
  output {
    Int A_out = read_int(stdout())
  }
 runtime { continueOnReturnCode: false }
}

task B {
  Int B_in
  command {
    echo ${B_in}
  }
}


workflow w {
  call A
  call B {input: B_in = A.A_out}
}
