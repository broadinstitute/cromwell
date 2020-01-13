
task foo {
  command {
    echo "foo"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf {
  scatter (i in range(2)) {
    call foo
  }
}
