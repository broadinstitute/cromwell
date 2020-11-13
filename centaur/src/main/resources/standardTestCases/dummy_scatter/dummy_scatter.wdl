version 1.0

workflow dummy_scatter {
  scatter (x in range(35000)) {
    call dummy_scattered_task
  }
  output {
    Int results_count = length(dummy_scattered_task.string_out)
  }
}

task dummy_scattered_task {
  command {
    echo hello
  }
  output {
    String string_out = "hello"
  }
  runtime {
    # This is technically unnecessary given the test setup, but I want to double check this isn't accidentally sent
    # to a _real_, money-spending backend:
    backend: "Dummy"
  }
}
