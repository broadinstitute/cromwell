task makeGcsFileInJes {
  command {
    echo "I am writing stuff" > out
    sleep 2
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    File outputFile = "out"
  }
}

task readGcsFileFromLocal {
  File gcsFile

  command {
    cat ${gcsFile}
    sleep 2
  }

  runtime {
    backend: "Local"
  }

  output {
    String result = read_string(stdout())
  }
}

workflow local_gcs {
  call makeGcsFileInJes
  call readGcsFileFromLocal { input: gcsFile = makeGcsFileInJes.outputFile }
}
