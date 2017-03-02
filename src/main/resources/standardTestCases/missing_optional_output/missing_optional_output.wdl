task x {
    Int i
    command { echo $i }
    runtime { docker: "ubuntu" }
    output {
        Int out_but_intentionally_misnamed = i
        Boolean validOutput = i % 2 == 0
    }
}

workflow missing_optional_output {
  Array[Int] arr = [0,1,2,3]
  scatter (i in arr) {
    call x { input: i = i }
    if (x.validOutput) {
      Int x_out = x.out_except_undeclared
    }
  }

  Array[Int?] x_out_maybes = x_out
}
