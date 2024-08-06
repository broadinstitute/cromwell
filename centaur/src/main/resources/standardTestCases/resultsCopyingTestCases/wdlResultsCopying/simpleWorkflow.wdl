workflow simpleWorkflow {
    call simpleStdoutTask {}
    output {
        File outFile = simpleStdoutTask.outFile
        Array[File] outGlob = simpleStdoutTask.outGlob
    }
}

task simpleStdoutTask {
  String outputFileName = "output.txt"

  command {
    echo 'Hello world' > ${outputFileName}
    echo 'foo' > "foo.zardoz"
    echo 'bar' > "bar.zardoz"
    echo 'baz' > "baz.zardoz"
  }

  runtime {
    docker: "ubuntu"
    memory: "2G"
    cpu: 1
  }

  output {
    File outFile = outputFileName
    Array[File] outGlob = glob("*.zardoz")
  }
}
