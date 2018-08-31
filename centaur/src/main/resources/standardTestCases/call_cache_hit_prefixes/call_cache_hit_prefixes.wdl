version 1.0

task yo {
  input {
    String salutation
  }
  command {
    echo 'sup ~{salutation}?'
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    String out = read_string(stdout())
  }
}

workflow call_cache_hit_prefixes {
  call yo

  output {
    String sup = yo.out
  }
}
