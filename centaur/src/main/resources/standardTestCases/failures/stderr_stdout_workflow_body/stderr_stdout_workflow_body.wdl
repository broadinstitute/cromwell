version 1.0

workflow break_with_stderr {

  output {
    File load_data_csv = select_first([stdout(), stderr()])
  }
}
