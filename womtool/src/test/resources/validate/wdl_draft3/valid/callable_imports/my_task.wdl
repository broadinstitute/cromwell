version draft-3

task my_task {
  input {
    String s
  }

  command <<<
    echo -n "~{s}" | wc | awk ' { print $3  } '
  >>>

  output {
    Int i = read_int(stdout())
  }
}
