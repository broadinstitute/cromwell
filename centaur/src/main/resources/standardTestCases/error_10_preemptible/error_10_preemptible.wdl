version 1.0

task delete_self_if_preemptible {

  command {

    preemptible=$(curl -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible")

    # Delete self if running on a preemptible VM. This should produce an "error 10" which Cromwell should treat as a preemption.
    # Since `preemptible: 1` the job should be restarted on a non-preemptible VM.
    if [ "$preemptible" = "TRUE" ]; then
    
      fully_qualified_zone=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
      zone=$(basename "$fully_qualified_zone")

      gcloud compute instances delete $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=$zone -q
    fi

  }

  runtime {
    preemptible: 1
    docker: "google/cloud-sdk:slim"
  }
}


workflow error_10_preemptible {
  call delete_self_if_preemptible
}
