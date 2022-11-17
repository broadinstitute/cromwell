version 1.0

task imitate_oom_error {
  command {
    printf "Exception in thread "main" java.lang.OutOfMemoryError: testing\n\tat Test.main(Test.java:1)\n" >&2
    touch foo
    exit 137
  }
  output {
    File foo = "foo"
  }
  runtime {
    docker: "python:latest"
    memory: "1 GB"
    maxRetries: 2
    backend: "Papiv2"
  }
}

workflow retry_with_more_memory_after_137 {
  call imitate_oom_error
}
