version 1.0

workflow sleepy_sleep {

  input {
    Int sleep_seconds = 180
  }

  call sleep {
    input: sleep_seconds = sleep_seconds
  }

}

task sleep {

  input {
    Int sleep_seconds
  }

  meta {
    volatile: true
  }

  # `northamerica-northeast1` is Montreal, which is is geographically close (low latency),
  # uses 100% carbon-free energy, and best of all... we don't use for anything else,
  # so I can drop its quota to something ridiculousy low in `broad-dsde-cromwell-dev`
  runtime {
    cpu: 12
    docker: "ubuntu:latest"
    zones: "northamerica-northeast1-a northamerica-northeast1-b northamerica-northeast1-c"
  }

  command <<<
    sleep ~{sleep_seconds};
    ls -la
  >>>
}
