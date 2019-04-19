version 1.0

workflow cycle {
  call mirror as m1 {
    input: i = m3.out
  }
  call mirror as m2 {
    input: i = m1.out
  }

  call mirror as m3 {
    input: i = m2.out
  }

  call mirror as m4 {
    input: i = m3.out
  }
}

task mirror {
  input {
    Int i
  }

  command <<<>>>
  output {
    Int out = i
  }
}
