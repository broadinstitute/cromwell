task print {
  command {
    echo "She sells sea shells by the sea shore"
  }
  output {
    File tongueTwister = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task grep {
  String match = "o"
  File input_file
  command {
    grep '${match}' ${input_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow writeToCache {
  call print
  call grep { input: input_file = print.tongueTwister }
  call print as printAgain
  call grep as grepAgain { input: input_file = print.tongueTwister }
  output {
    grep.count
    grepAgain.count
  }
}
