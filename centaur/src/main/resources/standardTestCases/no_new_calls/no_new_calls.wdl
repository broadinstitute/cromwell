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
    command {
     sleep 100
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
