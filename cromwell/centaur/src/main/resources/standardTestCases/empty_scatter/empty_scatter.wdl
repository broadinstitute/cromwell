workflow empty_scatter {
  Array[Int] xs = []

  scatter (x in xs) {
    Int decl = 75 + x
    call do_nothing { input: whoa = decl }
  }

  output {
    Int task_outs_length = length(do_nothing.str)
    Int decls_length = length(decl)
  }
}

# We should never invoke this since we're scattering over an empty array
task do_nothing {
  Int whoa
  command {
    # doesn't matter, won't be run anyway:
    exit 500
  }
  output {
    String str = "str"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
