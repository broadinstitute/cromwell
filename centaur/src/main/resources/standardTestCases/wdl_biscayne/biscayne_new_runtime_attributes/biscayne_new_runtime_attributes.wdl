version development-1.1

workflow runtime_attributes_wf {
  call runtime_attributes_task
  output {
    String out = runtime_attributes_task.out
  }
}

task runtime_attributes_task {

  command <<<
    echo "Zardoz"
  >>>

  runtime {
    # Meaningless keys are ignored
    "banana": object {
        cpuPlatform: "Banana Lake"
    }

    "gcp" : object {
        # Platform-specific keys take precedence
        docker: "rockylinux:9"
    }

    # Generic keys are ignored in favor of platform ones
    docker: "ubuntu:latest"

    # We still read generic keys that are not overridden
    cpu: 4
  }

  output {
    String out = read_string(stdout())
  }
}
