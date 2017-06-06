task noop {
  command {}
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow workflow_type_and_version_wdl {
  call noop
}
