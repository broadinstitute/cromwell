version development

workflow afters {
  input {
    String where = "/tmp/helloFile"
  }

  # Should not impact 'read' because the second write overwrites it:
  call write_to_shared as foo1 { input: i = 5, where = where }

  # Scatter once to replace the contents with something else:
  scatter (a in range(1)) {
    call write_to_shared as foo2 after foo1 { input: i = 6, where = where }
  }

  # The call to 'read':
  call read_from_shared after foo2 { input: where = where }

  # Should not impact 'read' because it happens afterwards:
  call write_to_shared as foo3 after read_from_shared { input: i = 77, where = where }

  output {
    Int result = read_from_shared.read_result
  }
}

task write_to_shared {
  input {
    Int i
    String where
  }
  command <<<
    sleep 2
    echo "~{i}" > /tmp/helloFile
  >>>
  runtime {
    backend: "LocalNoDocker"
  }
  output {
    Boolean done = true
  }
}

task read_from_shared {
  input {
    String where
    Boolean ready = true
  }
  command <<<
    cat "~{where}"
  >>>
  output {
    Int read_result = read_int(stdout())
  }
  runtime {
    backend: "LocalNoDocker"
  }
}
