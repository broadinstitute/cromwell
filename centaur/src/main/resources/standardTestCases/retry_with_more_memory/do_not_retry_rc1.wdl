version 1.0

task imitate_oom_error {
  command {
    # This task mimics printing an exception that would cause it to retry with increased memory, but
    # because the return code is 1 and continueOnReturnCode is true, the task does not retry.
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

workflow do_not_retry_rc1 {
  call imitate_oom_error
}
