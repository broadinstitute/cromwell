task A {
  String i
  command {ps}
  output {String o = read_string(stdout())}
}

workflow w {
  Int i
  Array[String] arr

  call A

  if (i == 2) {
    call A as B
  }

  if (A.o == "foo") {
    call A as C
  }

  if (A.o == "bar") {
    scatter(x in arr) {
      call A as D {input: i=x}
    }
  }

  call A as E {input: i=C.o}
}
