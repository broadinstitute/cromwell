task firstTask {
  command {
    echo foo
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}

task secondTask {
  String relay
  command {
    echo ${relay}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow sub_decls {
  # A subworkflow declaration that requires an input.
  String uninitialized
  # A subworkflow declaration that's initialized but which we expect to be overridden.
  String wantsOverride = "don't want to see this"

  call firstTask
  # This subworkflow declaration is initialized and not overridden, so it shouldn't be
  # evaluated at subworkflow call time. If it is evaluated at subworkflow call time that
  # will fail since the task call on which this depends hasn't run yet.
  String firstTaskOut = firstTask.out
  call secondTask { input: relay = firstTaskOut }
  output {
    String out = secondTask.out
    String init = uninitialized
    String override = wantsOverride
  }
}
