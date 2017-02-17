task print {
  command {
    echo "She sells sea shells by the sea shore only on $(date)"
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
    File redirect = input_file
  }
  runtime {docker: "ubuntu:latest"}
}

task grepAgain {
  String match = "o"
  File input_file
  command {
    grep '${match}' ${input_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
    File redirect = input_file
  }
  runtime {docker: "ubuntu:latest"}
}

workflow writeToCache {
  call print
  call grep { input: input_file = print.tongueTwister }
  call grepAgain { input: input_file = grep.redirect }
  output {
    grep.count
    grepAgain.count
  }
}
