workflow Test {
  File wfIn

  call TestTask { input: taskIn=wfIn }
}

task TestTask {
  File taskIn
  # Array[File] taskArray

  command {
    cat ${taskIn} > "out.txt"
  }

  output {
    File taskOut = "out.txt"
  }

  runtime {
    foo: "bar"
    docker: "bar"
    dockerWorkingDir: "foo"
  }
}
