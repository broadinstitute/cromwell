import "sub_workflow_interactions_import.wdl" as counter

task hello {
  String addressee
  
  command {
    echo "Hello ${addressee}!" > hello
    wc -w < hello > count
    sleep 2
  }
  runtime {
      docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string("hello")
    Int count = read_int("count")
  }
}

task read {
    Array[File] files
    command {
        cat ${sep=' ' files}
        sleep 2
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        String read_out = read_string(stdout())
    }
}

workflow sub_workflow_interactions_scatter {
  Array[Int] arr = [0, 2]
  
  scatter(i in arr) {
    call hello as helloInScatter { input: addressee = i }
    call counter.countEvens as countEvensInScatter { input: max = helloInScatter.count + i } # sub workflow in scatter
    call hello as secondHelloInScatter { input: addressee = countEvensInScatter.someStringOutput }
  }
  
  call read { input: files = countEvensInScatter.evenFile }
  
  output {
    # new output syntax
    Array[String] out = helloInScatter.salutation
    Array[String] out3 = secondHelloInScatter.salutation
    String out6 = read.read_out
  }
}