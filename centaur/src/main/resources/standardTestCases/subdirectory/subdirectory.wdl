task subdirTask {
  command {
    mkdir subdir
    cd subdir
    echo "I'm in a subdirectory !" > subFile
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    File outputFile = "subdir/subFile"
  }
}

workflow subdirWorkflow {
  call subdirTask
}
