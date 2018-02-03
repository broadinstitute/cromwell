# The idea of this test is that this workflow will fail to instantiate, but we want to make sure labels are populated prior to that error

task hello_task {
  command {
    echo "foo" > foo.txt
  }

  runtime {
    docker: "python:2.7"
  }

  output {
    File foo = "foo.txt"
  }
}

task hello_task2 {
  command {
    echo "foo" > foo.txt
  }

  runtime {
    docker: "python:2.7"
  }

  output {
    File foo = "foo.txt"
  }
}

workflow hello {
  String foo

  call hello_task {
    input:
      asdf=foo
  }
  call hello_task2

  output {
    File foo = hello_task.foo
    File foo2 = hello_task2.foo
  }
}
