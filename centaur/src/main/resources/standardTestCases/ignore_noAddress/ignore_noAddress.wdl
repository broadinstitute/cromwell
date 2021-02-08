task task_with_noAddress {

    command { echo "hello, world" }

    output { String s = read_string(stdout()) }

    runtime {
      docker: "ubuntu:latest"
      # NB: This is ignored because of the config in the centaur configuration.
      # If we didn't have that, this task would run forever - or at least until PAPI times it out.
      noAddress: true
    }
}

workflow ignore_noAddress {
    call task_with_noAddress
}
