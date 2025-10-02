version development-1.1

task sayHi {
  input {
    String name
    String greeting
  }
  command <<<
    cat "~{greeting} ~{name}" > /tmp/helloFile
  >>>
}

workflow callInputPassthrough {
  String name = "friend"
  String greeting = "hello"

  # These two inputs should be handled the same way
  call sayHi {
    input:
      name,
      greeting = greeting
  }
}
