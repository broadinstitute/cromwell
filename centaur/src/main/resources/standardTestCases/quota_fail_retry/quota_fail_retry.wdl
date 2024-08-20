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

  # I set `broad-dsde-cromwell-dev` to have super low CPU quota in `us-west3` (Salt Lake City) for this test
  runtime {
    cpu: 10
    docker: "ubuntu:latest"
    zones: "us-west3-a us-west3-b us-west3-c"
  }

  command <<<
    sleep ~{sleep_seconds};
    ls -la
  >>>
}
