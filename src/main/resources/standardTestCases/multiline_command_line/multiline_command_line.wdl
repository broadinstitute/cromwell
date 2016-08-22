task blah {
  command <<<
    python <<CODE
    def a():
      return "a"
    def b():
      return "b"
    print('{}{}'.format(a(),b()))
    CODE
  >>>

  output {
    String ab = read_string(stdout())
  }

  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf {
  call blah
}
