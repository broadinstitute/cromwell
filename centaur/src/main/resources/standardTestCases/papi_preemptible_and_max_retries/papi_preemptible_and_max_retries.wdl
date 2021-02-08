version 1.0

task delete_self {

  command {
      fully_qualified_zone=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
      zone=$(basename "$fully_qualified_zone")

      gcloud compute instances delete $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=$zone -q
  }

  runtime {
    preemptible: 1
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:slim"
    maxRetries: 1
  }
}

workflow papi_preemptible_and_max_retries {
  call delete_self
}
