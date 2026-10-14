version 1.2

# Invalid: the deprecated runtime section cannot be combined with requirements or hints
task runtime_with_requirements {
  command <<<
    echo hello
  >>>

  runtime {
    docker: "ubuntu:latest"
  }

  requirements {
    cpu: 1
  }
}

workflow runtime_with_requirements_wf {
  call runtime_with_requirements
}
