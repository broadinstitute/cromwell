task hey {
  String what
  command {
    echo "Hey ${what}!"
  }
  output {
    String lyrics = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow refresh_token_failure {
  call hey
  output {
     hey.lyrics
  }
}