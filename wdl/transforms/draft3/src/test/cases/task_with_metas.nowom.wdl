version draft-3

task task_with_metas {
  input {
    Int a
    Int b
  }

  output {
    Int out = read_int(stdin())
  }

  command {
    echo $((${a} + ${b}))
  }

  meta {
    author: "John Doe",
    email: "john.doe@yahoo.com"
  }

  parameter_meta {
    a: "just an integer",
    b: "an important parameter"
  }
}
