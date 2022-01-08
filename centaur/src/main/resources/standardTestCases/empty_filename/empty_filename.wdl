version 1.0

workflow empty_filename {
  input {
    File file
  }

  call use_file { input: file = file }
}

task use_file {
  input {
    File file
  }
  command <<<
    cat ~{file}
  >>>
  output {
    String content = read_string(stdout())
  }
}
