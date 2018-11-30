version development

workflow after_circle {
  call foo as foo1 after foo2 { input: i = 5 }
  call foo as foo2 after foo1 { input: i = 6 }
}

task foo {
  input {
    Int i
  }
  command <<<
    cat "hello ~{i}" > /tmp/helloFile
  >>>
}
