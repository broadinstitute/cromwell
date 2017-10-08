task cgrep {
  String pattern
  File in_file

  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}
