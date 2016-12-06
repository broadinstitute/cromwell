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

workflow sub_workflow_interactions {
  call hello { input: addressee = "Sub Workflow World" }
  call counter.countEvens { input: max = hello.count } # Sub worklfow depends on previous task call
  call hello as secondHello { input: addressee = countEvens.someStringOutput } # Task call depends on previous sub workflow call
  
  output {
    # old output syntax
    hello.* 
    
    # new output syntax
    String out2 = secondHello.salutation
    Array[Int] out4 = read_lines(countEvens.evenFile)
  }
}