version development

workflow afters_and_scatters {

  meta {
    description: "Makes sure that we can process 'after' definitions between scatters, within the same scatter, and from within an outer scatter"
  }

  input {
    String where = "/tmp/helloFile"
  }

  scatter (a in range(3)) {
    # Should not impact 'read' because foo2 overwrites it:
    call write_to_shared as foo0 { input: i = 0, where = "~{where}~{a}" }
    call read_from_shared as read0 after foo0 { input: where = "~{where}~{a}" }
  }

  # Scatter once to replace the contents with something else:
  scatter (a in range(3)) {
    # Should not impact 'read' because foo2 overwrites it:
    call write_to_shared as foo1  { input: i = 1, where = "~{where}~{a}" }

    call write_to_shared as foo2 after foo1 { input: i = 2 + a, where = "~{where}~{a}" }

    # The call to 'read':
    call read_from_shared after foo2 { input: where = "~{where}~{a}" }
  }

  # Should not impact 'read' because it happens afterwards:
  call write_to_shared as foo3 after read_from_shared { input: i = 7, where = "~{where}0" }

  # A call to 'read' again:
  call read_from_shared as read2 after foo3 { input: where = "~{where}0" }

  output {
    Array[Int] results = read_from_shared.read_result
    Int read2_result = read2.read_result
  }
}

task write_to_shared {
  input {
    Int i
    String where
  }
  command <<<
    sleep ~{i}
    echo "~{i}" > "~{where}"
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
