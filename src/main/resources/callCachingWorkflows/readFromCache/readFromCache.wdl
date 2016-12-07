task echo {
  command {
    echo "Peter Piper picked a peck of pickled peppers"
    sleep 2
  }
  output {
    File out = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task find {
  String match = "r"
  File in_file
  command {
    grep '${match}' ${in_file} | wc -l
    sleep 2
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow readFromCache {
  call echo
  call find { input: in_file = echo.out }
  call echo as echoAgain
  call find as findAgain { input: in_file = echo.out }
  output {
    find.count
    findAgain.count
  }
}
