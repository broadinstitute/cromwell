version development

workflow no_input_no_output {
  call no_inputs

  call no_inputs{}

  call no_inputs {}

  call no_inputs {

  }

  call no_inputs
  {
  }

  call no_inputs {
    # A comment does not mess up the works
  }

}

task no_inputs {
  command { echo Hello World }

}
