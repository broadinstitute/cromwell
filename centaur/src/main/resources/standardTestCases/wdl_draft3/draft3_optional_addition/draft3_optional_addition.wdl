version 1.0

workflow draft3_optional_addition {
  scatter(x in range(2)) {
    if (x % 2 == 0) {
      Int i = x
    }

    call call_me { input: maybe_i = i }
  }

  output {
    Array[String] lines = call_me.s
  }
}

task call_me {
  input {
    Int? maybe_i
  }
  command {
    echo "cmd~{" --only-with-maybe_i " + maybe_i}"
  }

  output {
    String s = read_string(stdout())
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
