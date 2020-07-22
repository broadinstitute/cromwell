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
}

workflow echo_wf {
  call string_interpolation_task
}
