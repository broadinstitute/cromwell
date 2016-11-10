workflow Test {
  File wfIn

  call TestTask { input: taskIn=wfIn }
}

task TestTask {
  File taskIn

  command {
    cat ${taskIn}
  }

  output {
    File taskOut = "out.txt"
  }

  runtime {
    docker: "bar"
  }
}
