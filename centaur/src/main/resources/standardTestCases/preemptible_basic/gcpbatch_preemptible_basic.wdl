version 1.0

task delete_self_if_preemptible {

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    preemptible=$(curl -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/scheduling/preemptible")

    # Perform a maintenance event on this VM if it is preemptible, which should cause it to be preempted.
    # Since `preemptible: 1` the job should be restarted on a non-preemptible VM.
    if [ "$preemptible" = "TRUE" ]; then
      fully_qualified_zone=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
      zone=$(basename "$fully_qualified_zone")

      gcloud beta compute instances simulate-maintenance-event $(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google") --zone=$zone -q
      sleep 60
    fi

  >>>

  runtime {
    preemptible: 1
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:slim"
  }
}


workflow preemptible_basic {
  call delete_self_if_preemptible
}
