task sum {
  File gumby
  command {
    /usr/bin/md5sum ${gumby} > file.md5
  }
  output {
    File out = "file.md5"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task cromwell_killer {
  File sequential_input
  command {
    sleep 60
    echo restarting yo
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow papi_upgrade {
  call sum
  call cromwell_killer { input: sequential_input = sum.out }
}
