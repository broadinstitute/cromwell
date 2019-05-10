version 1.0

workflow gpu_on_papi {

  input {
    Array[String] gpus
  }

  scatter (gpu in gpus) {
    call task_with_gpu { input: gpuTypeInput = gpu }
  }

  output {
    Array[Int] reported_gpu_counts = task_with_gpu.gpuCount
    Array[String] reported_gpu_types = task_with_gpu.gpuType
  }
}


task task_with_gpu {
  input {
    String gpuTypeInput
  }

  command {
    curl "https://www.googleapis.com/compute/v1/projects/broad-dsde-cromwell-dev/zones/us-central1-c/instances/$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")?fields=guestAccelerators" -H "Authorization: Bearer $(gcloud auth application-default print-access-token)" -H 'Accept: application/json' --compressed
  }

  output {
    Object metadata = read_json(stdout())
    Int gpuCount = metadata.guestAccelerators[0].acceleratorCount
    String gpuType = metadata.guestAccelerators[0].acceleratorType
  }

  runtime {
    gpuCount: 1
    gpuType: gpuTypeInput
    docker: "google/cloud-sdk:slim"
    zones: ["us-central1-c"]
  }
}
