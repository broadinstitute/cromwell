version 1.0

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
    x: ["A", "B", "C"]
    y: [1, 2, 3]
    yf: [1.1, 2.9, 3.14]
    z: {k1: 1, k2: 2, k3: 3}
  }
}
