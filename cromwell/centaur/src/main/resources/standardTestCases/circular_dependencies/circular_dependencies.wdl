task mirror {
  String a
  String b
  command {
    echo ${a}_${b}
  }
  output {
    String out = read_string(stdout())
  }
}

workflow circular_dependencies {
  String hello = "hello"

  call mirror as mirror1 { input: a = hello, b = mirror2.out }
  call mirror as mirror2 { input: a = hello, b = mirror1.out }
}
