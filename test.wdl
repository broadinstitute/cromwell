workflow Test {
  File wfIn

  call TestTask { input: taskIn=wfIn }
}

task TestTask {
  File taskIn

  command {
    cat ${taskIn}
  }

  runtime {
    docker: "bar"
  }
}
