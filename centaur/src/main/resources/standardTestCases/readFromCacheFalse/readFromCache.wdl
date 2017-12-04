task echo {
  command {
    echo "Peter Piper picked a peck of pickled peppers"
  }
  output {
    File out = stdout()
  }
  runtime {docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"}
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
  runtime {docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"}
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
