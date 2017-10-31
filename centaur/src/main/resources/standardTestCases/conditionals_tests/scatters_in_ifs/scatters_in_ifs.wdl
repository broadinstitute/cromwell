task mirror {
  Int i
  command {
    echo ${i}
  }
  output {
    Int out = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}


workflow scatters_in_ifs {

# TODO: Reinstate the nested lookup with #2724
#  Array[Int] numbers = range(3)

  if (true) {
    scatter (n in range(3)) {
      call mirror as mirrorTrue { input: i = n }
    }
  }

  if (false) {
    scatter (n in range(3)) {
      call mirror as mirrorFalse { input: i = n }
    }
  }

  output {
    Array[Int]? inTruth = mirrorTrue.out
    Array[Int]? inLies = mirrorFalse.out
  }
}
