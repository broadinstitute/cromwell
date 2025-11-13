version 1.0

task imitate_oom_error_on_preemptible {
  command {
    # This task simulates a scenario where the tool is running and almost running OOM but before it does, the VM is
    # preempted. On the next attempt, Cromwell will use a non-preemptible VM because `preemptible: 1` in runtime
    # attributes. Since the VM was preempted on first attempt, Cromwell should not increase the memory on for second attempt.
    # However, on second attempt, the VM won't preempt leading to a OOM exception simulation. And hence on the
    # third attempt, Cromwell should have retried with more memory.

    echo "Doing meaningful work..."

    preemptible=$(curl -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible")

    # Simulate a maintenance event on ourselves if running on a preemptible VM
    # Since `preemptible: 1` the job should be restarted on a non-preemptible VM.
    if [ "$preemptible" = "TRUE" ]; then
      fully_qualified_zone=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
      zone=$(basename "$fully_qualified_zone")

      gcloud beta compute instances simulate-maintenance-event $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=$zone -q
      sleep 300
    fi

    # Should reach here on the second attempt
    printf "Exception in thread "main" java.lang.OutOfMemoryError: testing\n\tat Test.main(Test.java:1)\n" >&2
    exit 1
  }

  runtime {
    preemptible: 1
    docker: "mirror.gcr.io/google/cloud-sdk:slim"
    memory: "1 GB"
    maxRetries: 1
    backend: "GCPBATCH"
  }
}

workflow preemptible_and_memory_retry {
  call imitate_oom_error_on_preemptible
}
