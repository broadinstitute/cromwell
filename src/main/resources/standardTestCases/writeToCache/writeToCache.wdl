task print {
  command {
    echo "She sells sea shells by the sea shore only on $(date)"
  }
  output {
    File tongueTwister = stdout()
  }
  runtime {docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"}
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
  runtime {docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"}
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
  runtime {docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"}
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
