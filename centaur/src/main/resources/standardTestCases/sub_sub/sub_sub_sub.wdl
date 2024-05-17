
task foo {
  command {
    echo "foo"
    exit 1
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
