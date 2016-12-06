task countTo {
    Int value
    command {
        seq 0 1 ${value}
        sleep 2
    }
    runtime {
          docker: "ubuntu:latest"
      }
    output {
        File range = stdout()
    }
}

task filterEvens {
    File numbers
    command {
        grep '[02468]$' ${numbers} > evens
        sleep 2
    }
    runtime {
          docker: "ubuntu:latest"
      }
    output {
        File evens = "evens"
    }
}

workflow countEvens {
    Int max = 10
    
    call countTo { input: value = max }
    call filterEvens { input: numbers = countTo.range }
    output {
        String someStringOutput = "I'm an output"
        File evenFile = filterEvens.evens
    }
}