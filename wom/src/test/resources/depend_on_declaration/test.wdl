task a {
  command { ... }
  output { Int o = read_int(stdout()) }
}

task b {
  Array[Int] ints
  command { ... }
  output { Int o = read_int(stdout()) }
}

workflow w {
  call a as a0
  call a as a1
  Array[Int] ints = [a0.o, a1.o]
  call b {input: ints=ints}
}
