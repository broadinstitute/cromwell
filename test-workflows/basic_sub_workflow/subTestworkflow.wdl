version 1.0

task subSubworkflowHello {
  command {
    echo "${2/2}"
  }
  output {
    String out = read_string(stdout())
  }
}


task subSubworkflowGoodbye {
  command {
    echo "${2/0}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow sub_wf_scattering {
  Array[Int] numbers = [1, 2, 3]
  call subSubworkflowHello
  scatter (i in numbers) {
    call subSubworkflowGoodbye
  }
  output {
    String out = subSubworkflowHello.out
    Array[String] out2 = subSubworkflowGoodbye.out
  }
}
