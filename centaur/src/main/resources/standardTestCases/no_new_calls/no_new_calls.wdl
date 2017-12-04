task boundToFail {
  command {
    cd fake/file/path
  }
  output {
   String badOutput = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task shouldNotStart {
  String str
    command {
     echo ${str}
    }
    runtime {
       docker: "ubuntu:latest"
    }
}

task shouldSucceed {
  String str

  # For this test to pass, currently `sleep` must be large as JES can take minutes to run a job. Otherwise, if
  # shouldSuceed finishes before boundToFail, delayedTask2 may still attempt a restart.

  # Example where delayedTask2 attempts restart by 17:13:43, while boundToFail doesn't return until 13 seconds later:

  # boundToFail         submitted       17:08:14  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L3312
  # shouldSucceed       submitted       17:08:26  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L3616
  # boundToFail        initializing     17:08:27  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L3680
  # shouldSucceed      initializing     17:08:38  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L3735
  # shouldSucceed        running        17:09:13  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L3864
  # <cromwell>           restart        17:10:59  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L4825
  # boundToFail          restart        17:11:23  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L5170
  # shouldSucceed        restart        17:11:24  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L5270
  # boundToFail          running        17:11:34  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L5471
  # shouldSucceed        finish         17:11:35  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L5490
  # delayedTask2   >>    restart    <<  17:13:43  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L6212
  # boundToFail          finish         17:13:56  https://travis-ci.org/broadinstitute/cromwell/jobs/310544199#L6259

    command {
     sleep 300
     echo ${str}
    }
    runtime {
       docker: "ubuntu:latest"
    }
    output {
     String stalling = read_string(stdout())
    }
}

task delayedTask {
  String str_2
    command {
     echo ${str_2}
    }
    runtime {
       docker: "ubuntu:latest"
    }
    output {
     String delayedOut = read_string(stdout())
    }
}

workflow no_new_calls {
  call boundToFail
  call shouldNotStart { input: str = boundToFail.badOutput }
  call shouldSucceed { input: str = "echoing nonsense" }
  call delayedTask as delayedTask1 { input: str_2 = shouldSucceed.stalling }
  call delayedTask as delayedTask2 { input: str_2 = delayedTask1.delayedOut }
}
