version 1.0

workflow google_labels {

  meta {
    description: "Confirms that the google_labels workflow option is propagated correctly to tasks"
  }

  input {
    Array[Pair[String, String]] expected_kvps = [ ("wdl-task-name", "check-labels") ]
    Boolean check_aliases = true
  }

  # Build an array of checker lines to confirm the existence of each KVP:
  scatter (p in expected_kvps) {
    String checker_lines = "echo $INSTANCE_INFO | grep -q '\"~{p.left}\": \"~{p.right}\"' || echo \"Couldn't find correct '~{p.left}' label in instance info: $INSTANCE_INFO\" 1>&2"
  }

  call check_labels { input: label_checker_lines = checker_lines }

  if (check_aliases) {
    call check_labels as check_label_alias { input: label_checker_lines = [
      "echo $INSTANCE_INFO | grep -q '\"wdl-task-name\": \"check-labels\"' || echo \"Couldn't find correct 'wdl-task-name' label in instance info: $INSTANCE_INFO\" 1>&2",
      "echo $INSTANCE_INFO | grep -q '\"wdl-call-alias\": \"check-label-alias\"' || echo \"Couldn't find correct 'wdl-call-alias' label in instance info: $INSTANCE_INFO\" 1>&2"
    ] }
  }

  output {
    Boolean done = true
  }

}

task check_labels {

  input {
    Array[String] label_checker_lines
  }

  command {
    # Check that the task metadata is correct:
    INSTANCE=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
    TOKEN=$(gcloud auth application-default print-access-token)
    INSTANCE_INFO=$(curl -s "https://www.googleapis.com/compute/v1/projects/broad-dsde-cromwell-dev/zones/us-central1-c/instances/$INSTANCE" \
      -H "Authorization: Bearer $TOKEN" \
      -H 'Accept: application/json')

    # Check for the custom label in the instance info:
    ~{sep="\n" label_checker_lines}
  }
  runtime {
    docker: "google/cloud-sdk:slim"
    # Specify this zone because it's used in the curl commands above. We probably *could* work this out ad-hoc but it's easier to hard-code it here:
    zones: ["us-central1-c"]
    failOnStderr: true
  }
}
