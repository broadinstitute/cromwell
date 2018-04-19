workflow retry_for_me {

  call broken_task {}
}

task broken_task {

    String s = "/not/a/file"

    command { cat ${s} }
    runtime {
       docker: "ubuntu"
       maxRetries: 1
    }
}
