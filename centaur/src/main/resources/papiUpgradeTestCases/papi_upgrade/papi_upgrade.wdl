task sum {
  File gumby
  command {
    /usr/bin/md5sum ${gumby} > file.md5
  }
  output {
    File out = "file.md5"
  }
  runtime {
    docker: "ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90"
  }
}

task cromwell_killer {
  File sequential_input
  command {
    sleep 60
    echo restarting yo
  }
  runtime {
    docker: "ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90"
  }
}

workflow papi_upgrade {
  call sum
  call cromwell_killer { input: sequential_input = sum.out }
}
