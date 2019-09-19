version 1.0

task imitate_oom_error {
  command {
    printf "Exception in thread "main" java.lang.OutOfMemoryError: testing\n\tat Test.main(Test.java:1)\n" >&2 && (exit 1)
  }
  runtime {
    docker: "python:latest"
    memory: "1 GB"
    continueOnReturnCode: true
    maxRetries: 1
    backend: "Papiv2-Retry-With-More-Memory"
  }
}

workflow retry_with_more_memory {
  call imitate_oom_error
}
