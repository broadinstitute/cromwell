version development

workflow default_default {
  call default_default_task
  output {
    String out = default_default_task.out
  }
}

task default_default_task {
  input {
    String? x
    Int? y
  }

  command <<<
    echo "x: '~{x}' y: '~{y}'"
  >>>
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}
