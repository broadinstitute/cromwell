version 1.0

task delete_self {

  command {
      preemptible=$(curl -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible")

      # Simulate a maintenance event on ourselves if running on a preemptible VM, otherwise delete ourselves.
      fully_qualified_zone=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
      zone=$(basename "$fully_qualified_zone")

      if [ "$preemptible" = "TRUE" ]; then
        gcloud beta compute instances simulate-maintenance-event $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=$zone -q
      else
        # We need to actually delete ourselves if the VM is not preemptible; simulated maintenance events don't seem to
        # precipitate the demise of on-demand VMs.
        gcloud compute instances delete $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=$zone -q
      fi
      sleep 120
  }

  runtime {
    preemptible: 1
    docker: "mirror.gcr.io/google/cloud-sdk:slim"
    maxRetries: 1
  }
}

workflow papi_preemptible_and_max_retries {
  call delete_self
}
