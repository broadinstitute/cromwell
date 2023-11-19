version 1.0

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
  input {
    Int number
  }
  command {
    echo "${number/0}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow wf_scattering {
  Array[Int] numbers = [1, 2, 3]
  Array[Int] numbers2 = [1, 2]
  scatter (i in numbers) {
     call scatteringHello {
      input: number = i
     }
  }
  scatter (i in numbers2) {
    call scatteringGoodbye {
      input: number = i
    }
  }
  output {
    Array[String] out = scatteringHello.out
    Array[String] out2 = scatteringGoodbye.out
  }
}

