task blah {
  command <<<
    python <<CODE
    def a():
      return "a"
    def b():
      return "b"
    print('{}{}'.format(a(),b()))
    CODE
    sleep 2
  >>>

  output {
    String ab = read_string(stdout())
  }

  runtime {
    docker: "python:latest"
  }
}

workflow wf {
  call blah
}
