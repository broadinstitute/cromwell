version 1.0

import "subTestworkflow.wdl" as subTest

task scatteringHello {
  input {
    Int number
  }
  command {
    echo "${100/number}"
  }
  output {
    String out = read_string(stdout())
  }
}

task scatteringGoodbye {
  command {
    echo "${2/0}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow wf_scattering {
  Array[Int] numbers = [1, 2, 3]
  Array[Int] numbers2 = [1, 2]
  scatter (i in numbers) {
     call scatteringHello
  }
  scatter (i in numbers2) {
     call scatteringGoodbye {
      input: 
        number = 0
     }
  }
  call subTest.sub_wf_scattering
  output {
    Array[String] out = scatteringHello.out
    Array[String] out2 = scatteringGoodbye.out
  }
}

