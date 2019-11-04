
task foo {
  command {
    echo "foo"
  }
}

workflow wf {
  scatter (i in range(2)) {
    call foo
  }
}
