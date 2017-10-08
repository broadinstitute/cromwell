String global = "global"

task taskA {
  String i
  command {sh script.sh ${i}}
  output {String o = read_string(stdout())}
}

task taskB {
  Int i
  command {python add_one.py ${i}}
  output {Int o = read_int(stdout())}
}

workflow w {
  Int i = 2
  String foo = "foo"
  Array[String] arr

  call taskA as A

  if (i == 2) {
    Int j = i + 1
    Int h = 100

    call taskB as add_first {input: i=j}
    call taskB as add_second {input: i=i+h}
    call taskA as B

    scatter (y in arr) {
      String k = y + foo

      if (y == "foo") {
        Int m = i + 200

        call taskB as add_third {input: i=m}
        call taskA as C {input: i=k+global}
      }
    }
  }

  call taskA as D {input: i=C.o[0]}
}
