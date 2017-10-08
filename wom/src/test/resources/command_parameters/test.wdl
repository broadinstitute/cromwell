task param_test {
  String a
  String b
  Array[String] c
  Int? d
  Array[Int]+ e
  Boolean f

  command <<<
  ./binary ${a} ${"-p " + b} ${sep="," c} ${default=9 d} ${sep="\t" e} ${true="--true" false="--false" f}
  >>>
}

workflow wf {
  call param_test
}
