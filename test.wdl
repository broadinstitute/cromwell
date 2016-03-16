task t {
  String s
  command {
    echo "stdout ${s}"
    echo "stderr" 1>&2
  }
  output {
    String err = read_string(stderr())
    String out = read_string(stdout())
  }
}

workflow w {
  call t
}
