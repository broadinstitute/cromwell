version 1.2

task task_variable {
  command <<<
    echo "~{task.name}" > name.txt
  >>>

  output {
    String name = read_string("name.txt")
    String description = task.meta.description
    String? container = task.container
  }

  requirements {
    container: "ubuntu:latest"
  }

  meta {
    description: "Reads the task variable"
  }
}

workflow task_variable_wf {
  call task_variable

  output {
    String name = task_variable.name
    String description = task_variable.description
    String? container = task_variable.container
  }
}
