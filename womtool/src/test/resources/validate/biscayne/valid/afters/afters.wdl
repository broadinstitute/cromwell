version development

workflow afters {
  call foo { input: i = 5 }
  call foo as foo2 after foo { input: i = 6 }
}

task foo {
  input {
    Int i
  }
  command <<<
    sleep 10
    echo "hello ~{i}" > /tmp/helloFile
  >>>
}
