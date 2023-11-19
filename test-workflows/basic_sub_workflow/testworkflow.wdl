version 1.0

task scatteringHello {
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
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
  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }
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

