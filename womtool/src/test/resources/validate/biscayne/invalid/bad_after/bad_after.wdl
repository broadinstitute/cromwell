version development-1.1

workflow bad_after {
  call foo { input: i = 5 }
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
