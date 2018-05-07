version 1.0

workflow no_input_no_output {
  call no_inputs
}

task no_inputs {
  command { echo Hello World }

}
