version 1.0

task imitate_oom_error {
  meta {
    volatile: true
  }
  command {
    echo "$MEM_SIZE $MEM_UNIT"

    # Current bashes do not do floating point arithmetic, Python to the rescue.
    LESS=$(python -c "print($MEM_SIZE < 1.21)")

    if [[ "$LESS" = "True" ]]
    then
      printf "Exception in thread "main" java.lang.OutOfMemoryError: testing\n\tat Test.main(Test.java:1)\n" >&2
      exit 1
    fi

    echo "$MEM_SIZE $MEM_UNIT" > memory_output.txt
  }
  output {
    String memory_output = read_string("memory_output.txt")
  }
  runtime {
    docker: "python:latest"
    memory: "1 GB"
    maxRetries: 2
    backend: "GCPBATCH"
  }
}

workflow retry_with_more_memory {
  call imitate_oom_error

  output {
    String memory_output = imitate_oom_error.memory_output
  }
}
