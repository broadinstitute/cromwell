version 1.0

task imitate_oom_error {
  command {
    printf "Exception in thread "main" java.lang.OutOfMemoryError: testing\n\tat Test.main(Test.java:1)\n" >&2 && (exit 1)
  }
  output {
  }
  runtime {
    docker: "python:latest"
    memory: "1 GB"
    continueOnReturnCode: true
    maxRetries: 2
    backend: "Papiv2"
  }
}

workflow do_not_retry_with_more_memory_rc1 {
  call imitate_oom_error
}
