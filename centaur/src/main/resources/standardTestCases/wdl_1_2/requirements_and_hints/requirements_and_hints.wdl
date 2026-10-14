version 1.2

task requirements_and_hints {
  input {
    File input_file
    Int cpus = 1
  }

  command <<<
    # The input file is localization optional, so only its path is used here
    echo "~{input_file}" > path.txt
    # Exit with a return code that the requirements section accepts
    exit 3
  >>>

  output {
    String path = read_string("path.txt")
  }

  requirements {
    container: "ubuntu:latest"
    cpu: cpus
    memory: "1 GiB"
    max_retries: 0
    return_codes: [0, 3]
  }

  hints {
    short_task: true
    inputs: input {
      input_file: hints {
        localization_optional: true
      }
    }
  }
}

workflow requirements_and_hints_wf {
  input {
    File input_file
  }

  call requirements_and_hints { input_file = input_file }

  output {
    String path = requirements_and_hints.path
  }
}
