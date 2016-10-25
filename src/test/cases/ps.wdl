task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}
