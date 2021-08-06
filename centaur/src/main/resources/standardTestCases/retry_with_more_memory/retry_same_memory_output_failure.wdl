version 1.0

task imitate_oom_error {
  command {
    # This task mimics printing an exception that would cause it to retry with increased memory, but
    # because the return code is 1 and continueOnReturnCode is true, the task does not retry with more memory.
    printf "Exception in thread "main" java.lang.OutOfMemoryError: testing\n\tat Test.main(Test.java:1)\n" >&2 && (exit 1)
    # As a simulation of an OOM condition, do not create the 'foo' file. Cromwell should still be able to delocalize important detritus.
    # touch foo
  }
  output {
    # Since the file does not exist, this will cause the status to be Failed. The task will end up retrying,
    # but memory does not increase.
    File foo = "foo"
  }
  runtime {
    docker: "python:latest"
    memory: "1 GB"
    continueOnReturnCode: true
    maxRetries: 2
    backend: "Papiv2"
  }
}

workflow retry_same_memory_output_failure {
  call imitate_oom_error
}
