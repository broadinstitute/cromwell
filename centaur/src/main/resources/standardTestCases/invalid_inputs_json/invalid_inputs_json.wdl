task if_i_had_valid_json_inputs {
  command {
    echo "I would be so happy"
  }
}

workflow invalid_input_json {
  call if_i_had_valid_json_inputs
}
