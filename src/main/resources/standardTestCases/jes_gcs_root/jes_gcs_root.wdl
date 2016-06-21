#The associated workflow option has a specific gs path as the jes root, and this task uses outputs to confirm root.
task checkDir {
    command {
        cd..
        pwd
    }
    output {
        String root = read_string(stdout())
    }
    runtime {
    docker: "ubuntu:latest"
    }
}
workflow jes_wdl {
  call checkDir
  output {
    checkdir.root
  }
}
