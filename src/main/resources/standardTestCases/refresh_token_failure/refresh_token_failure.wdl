task hey {
  String what
  command {
    echo "Hey ${what}!"
    sleep 2
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