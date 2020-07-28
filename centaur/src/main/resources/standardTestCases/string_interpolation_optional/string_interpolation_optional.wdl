version development
task string_interpolation_task {

  input {
    String? str
  }

  command <<<
      echo ~{'"' + str + '"'}
  >>>

  output {
    String out = read_string(stdout())
  }
  
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow echo_wf {
  call string_interpolation_task
}
