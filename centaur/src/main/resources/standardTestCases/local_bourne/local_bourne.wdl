task shell {
  command {
    # The SHELL variable isn't set in /bin/sh on alpine, busybox or ubuntu so there is this abomination instead.
    ps auxwww | head -n 2 | tail -n 1 | cut -c22-28
  }
  output {
    String out = read_string(stdout())
  }
  runtime {
    docker: "busybox:1.28.3"
    backend: "LocalBourneShell"
  }
}

workflow local_bourne {
  call shell
  output {
     shell.out
  }
}
