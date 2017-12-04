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
    docker: "ubuntu"
  }
}
