
version 1.0

workflow file_inputs_in_output_expressions {
  call mkFoo
  call foo2bar { input: foo = mkFoo.foo, foo_index = mkFoo.foo_index  }
  call check_bar { input: bar = foo2bar.bar, bar_index = foo2bar.bar_index }
}

task mkFoo {
  command <<<
    echo "foo" >> foo.txt
    echo "fooidx" >> foo.txt.index
  >>>

  output {
    File foo = "foo.txt"
    # File foo_index = "~{foo}.index"
    File foo_index = "foo.txt.index"

  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task foo2bar {
  meta {
    # volatile so that we never call cache this. Useful when debugging, irrelevant in CI situations.
    volatile: true
  }

  input {
    File foo
    File foo_index
  }

  command <<<
    cp ~{foo} ~{"~{foo}.barre"}
    echo "bar" >> ~{foo}.barre
    cp ~{foo_index} ~{"~{foo_index + ".barre"}"}
    echo "baridx" >> ~{foo_index}.barre
  >>>

  output {
    File bar = "~{foo}.barre"
    File bar_index = foo_index + ".barre"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task check_bar {
  meta {
    # volatile so that we never call cache this. Useful when debugging, irrelevant in CI situations.
    volatile: true
  }

  input {
    File bar
    File bar_index
  }

  command <<<
    set -e

    grep "foo" ~{bar}
    grep "bar" ~{bar}
    grep "fooidx" ~{bar_index}
    grep "baridx" ~{bar_index}
  >>>

  output {
    # No outputs. If the check passes then the task will succeed.
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
