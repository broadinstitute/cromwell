import "sub_workflow_interactions_import.wdl" as counter

task hello {
  String addressee
  
  command {
    echo "Hello ${addressee}!" > hello
    wc -w < hello > count
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
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        String read_out = read_string(stdout())
    }
}

workflow sub_workflow_interactions {
  call hello { input: addressee = "Sub Workflow World" }
  call counter.countEvens { input: max = hello.count } # Sub worklfow depends on previous task call
  call hello as secondHello { input: addressee = countEvens.someStringOutput } # Task call depends on previous sub workflow call
  
  scatter(i in read_lines(countEvens.evenFile)) {
    call hello as helloInScatter { input: addressee = i }
    call counter.countEvens as countEvensInScatter { input: max = helloInScatter.count + i } # sub workflow in scatter
    call hello as secondHelloInScatter { input: addressee = countEvensInScatter.someStringOutput }
  }
  
  call read { input: files = countEvensInScatter.evenFile }
  
  output {
    # old output syntax
    hello.* 
    
    # new output syntax
    Array[String] out = helloInScatter.salutation
    String out2 = secondHello.salutation
    Array[String] out3 = secondHelloInScatter.salutation
    Array[Int] out4 = read_lines(countEvens.evenFile)
    String out6 = read.read_out
  }
}