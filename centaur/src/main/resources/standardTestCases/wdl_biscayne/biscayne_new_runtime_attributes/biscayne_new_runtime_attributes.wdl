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

  meta {
    volatile: true
  }

  runtime {
    # Meaningless keys are ignored
    banana: object {
        cpuPlatform: "Banana Lake"
    }

    gcp: object {
        # Platform-specific keys take precedence
        docker: "rockylinux:9",
        memory: "6 GB"
    }

    azure: object {
        memory: "4 GB",
        docker: "debian:latest"
    }

    # Generic keys are ignored in favor of platform ones
    docker: "ubuntu:latest"
    memory: "8 GB"

    # We still read generic keys that are not overridden
    #
    # FIXME/TODO:
    #
    # Runtime attributes are casted to String somewhere between WDL parsing and metadata checking in Centaur.
    # As a temporary measure to make this test pass, this test is expecting the CPU value to be a String.
    # This is an issue in Centaur and does not affect production.
    cpu: 4
  }

  output {
    String out = read_string(stdout())
  }
}
