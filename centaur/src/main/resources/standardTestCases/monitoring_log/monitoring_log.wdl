workflow monitoring_log {
  call get_stats
}

task get_stats {
  command {
    sleep 50
  }
  output {
    Array[String] stats = read_lines("monitoring.log")
  }
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
}
