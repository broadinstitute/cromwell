version 1.0

workflow attempt_resilience {
  call maybe_preempted as definitely_preempted { input: preempt = true }
  call maybe_preempted as unlikely_preempted { input: preempt = false }

  # We're not really interested in the workflow outputs, we care more about accessing the call level metadata
  output {
    Int i = 1
  }
}

task maybe_preempted {
  input {
    Boolean preempt
  }

  meta {
    volatile: true
  }

  command <<<

    # Test whether we're on attempt 1 by checking whether this VM is preemptible:
    preemptible=$(curl -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible")

    # If this VM is preemptible, and we want to preempt then simulate preemption:
    if [[ "$preemptible" = "TRUE" && "~{preempt}" = "true" ]]; then
      INSTANCE_NAME=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
      fully_qualified_zone=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
      zone=$(basename "$fully_qualified_zone")
      gcloud compute instances delete ${INSTANCE_NAME} --zone=$zone -q
    fi

    echo "foo" > foo_out.txt
  >>>

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:slim"
    preemptible: 1
  }

  output {
    String foo_out = read_string("foo_out.txt")
  }
}
