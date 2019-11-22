version 1.0

task imitate_oom_error {
  command {
    printf "Exception in thread "main" java.lang.OutOfMemoryError: testing\n\tat Test.main(Test.java:1)\n" >&2 && (exit 1)
    # As a simulation of an OOM condition, do not create the 'foo' file. Cromwell should still be able to delocalize important detritus.
    # touch foo
  }
  output {
    File foo = "foo"
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
