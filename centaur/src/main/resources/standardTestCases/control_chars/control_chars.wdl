version 1.0

task control {
  command {
    echo "“Control characters should work with Carbonited metadata” — Cromwell"
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
}

workflow control_chars {
  call control
}
