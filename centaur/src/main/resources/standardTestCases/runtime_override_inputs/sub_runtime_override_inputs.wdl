version development-1.1

task foo {
  input {
    String? bar
  }

  String baz = select_first([bar, "baz"])

  command {
    echo "${baz}"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf {

  meta {
    allowNestedInputs: true
  }

  scatter (i in range(2)) {
    call foo
  }
}
