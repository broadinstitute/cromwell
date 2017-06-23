task size_task {
  Float sz

  command {
    echo "file has size ${sz}"
  }
  output {
    File out = stdout()
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow size_wf {
  File file
  call size_task { input: sz = size(file) }

  output {
  }
}
