version 1.0

task add_one {
  input {
    Int x
  }
  command <<<
    python -c "print({x} + 1)" > result.txt
  >>>
  output {
    Int y = read_int("result.txt")
  }
}

workflow simple {
  input {
    Int i
  }
  call add_one { input: x = i }
  output {
    Int result = add_one.y
  }
}