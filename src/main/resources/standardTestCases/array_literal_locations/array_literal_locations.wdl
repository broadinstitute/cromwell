# A task which takes an input, allows this to be called using an array literal in the input block
task array_literal_locations_i {
    Array[Int] array
    command { 
      echo ${sep=":" array} > out
      sleep 2
    }
    runtime { docker: "ubuntu:latest" }
    output { File out = "out" }
}

task array_literal_locations_ii {
    Int x
    command { 
      echo ${x} > out
      sleep 2
    }
    runtime {
        docker: "ubuntu:latest"
        continueOnReturnCode: [ 0, 1, 2 ]
    }
    output { File out = "out" }
}

workflow array_literal_locations {
    call array_literal_locations_i { input: array = [ 0, 1, 2 ] }
    scatter (x in [ 0, 1, 2 ]) {
        call array_literal_locations_ii { input: x = x }
    }
}
