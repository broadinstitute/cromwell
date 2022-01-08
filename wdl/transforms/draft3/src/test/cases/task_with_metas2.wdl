version 1.0

task task_with_metas2 {
  input {
  }

  output {
  }

  command {
  }

  meta {
    author: "John Doe"
    email: "john.doe@yahoo.com"
    b : true
    zipcode: 94043
    f : 1.3
    numbers : [1, 2, 3]
    extras: {
      house : "With porch",
      cat : "Lucy" }
  }
}
