version draft-3

task task_with_metas {
  input {
  }

  output {
  }

  command { echo Hello World }

  meta {
    author: "John Doe"
    email: "john.doe@yahoo.com"
  }

  parameter_meta {
    a: "just an integer"
    b: "an important parameter"
  }
}
