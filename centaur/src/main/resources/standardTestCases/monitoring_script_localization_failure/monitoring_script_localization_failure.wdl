task t {
  command {
    echo hi
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow monitoring_script_localization_failure {
  call t
}
