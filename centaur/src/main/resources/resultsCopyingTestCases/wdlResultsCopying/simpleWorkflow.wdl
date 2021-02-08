workflow simpleWorkflow {
    call simpleStdoutTask {}
    output {
        File outFile = simpleStdoutTask.outFile
    }
}

task simpleStdoutTask {
  String outputFileName = "output.txt"

  command {
    echo 'Hello world' > ${outputFileName}
  }

  runtime {
    docker: "ubuntu"
    memory: "2G"
    cpu: 1
  }

  output {
    File outFile = outputFileName
  }
}

