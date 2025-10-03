version development-1.1

task sayHi {
  input {
    String name
    String greeting
    String inquiry
  }
  command <<<
    cat "~{greeting} ~{name} ~{inquiry}" > /tmp/helloFile
  >>>
}

workflow callInputPassthrough {
  String name = "friend"
  String greeting = "hello"
  String inquiry = "how are you?"

  # These two inputs should be handled the same way
  call sayHi {
    input:
      name,
      greeting = greeting,
      inquiry
  }
}
