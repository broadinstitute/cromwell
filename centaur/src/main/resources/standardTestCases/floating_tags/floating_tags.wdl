task echo {
  command {
    echo "Peter Piper picked a peck of pickled peppers"
  }
  output {
    File out = stdout()
  }
  runtime {docker: "ubuntu:precise-20161209"}
}

task find {
  String match = "r"
  File in_file
  command {
    grep '${match}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:precise-20161209"}
}

workflow floating_tags {
  call echo
  call find { input: in_file = echo.out }
  call echo as echoAgain
  call find as findAgain { input: in_file = echo.out }
  output {
    find.count
    findAgain.count
  }
}
