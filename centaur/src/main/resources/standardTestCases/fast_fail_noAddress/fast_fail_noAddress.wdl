task task_with_noAddress {

    command { echo "hello, world" }

    output { String s = read_string(stdout()) }

    runtime {
      docker: "ubuntu:latest"
      # Our network is not configured specially, so this should cause this task to fail almost immediately.
      noAddress: true
    }
}

workflow fast_fail_noAddress {
    call task_with_noAddress
}
