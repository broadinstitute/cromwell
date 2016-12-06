import "refresh_token_sub_workflow_import.wdl" as sub

task hello {
  File iFile
  String addressee = read_string(iFile)
  command {
    echo "Hello ${addressee}!"
    sleep 2
  }
  output {
    String salutation = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow refresh_token_sub_workflow {
  call hello
  call sub.wf_hello { input: wfIFile = "gs://centaur-refresh-private/inputs/taylorSwift.txt" }
  output {
     hello.salutation
     wf_hello.salutation
  }
}
