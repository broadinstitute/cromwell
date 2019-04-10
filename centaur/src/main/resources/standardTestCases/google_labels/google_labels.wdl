version 1.0

workflow google_labels {
  call check_labels
}

task check_labels {

  meta {
    description: "Confirms that the google_labels workflow option is propogated correctly to tasks"
  }

  command {
    # Check that the task metadata is correct:
    INSTANCE=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
    TOKEN=$(gcloud auth application-default print-access-token)
    INSTANCE_INFO=$(curl -s "https://www.googleapis.com/compute/v1/projects/broad-dsde-cromwell-dev/zones/us-central1-c/instances/$INSTANCE" \
      -H "Authorization: Bearer $TOKEN" \
      -H 'Accept: application/json')

    # Check for the custom label in the instance info:
    echo $INSTANCE_INFO | grep -q '"custom-label": "custom-value"' || echo "Couldn't find custom-label in instance info: $INSTANCE_INFO" 1>&2
  }
  runtime {
    docker: "google/cloud-sdk:slim"
    zones: ["us-central1-c"]

    failOnStderr: true
  }
}
