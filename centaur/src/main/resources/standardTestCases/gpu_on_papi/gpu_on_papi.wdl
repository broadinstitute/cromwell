task task_with_gpu {
  String gpuTypeInput

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

workflow gpu_on_papi {
  call task_with_gpu as task_with_tesla_k80 { input: gpuTypeInput = "nvidia-tesla-k80" }
  call task_with_gpu as task_with_tesla_p100 { input: gpuTypeInput = "nvidia-tesla-p100" }

  output {
    Int tesla80GpuCount  = task_with_tesla_k80.gpuCount
    String tesla80GpuType = task_with_tesla_k80.gpuType
    Int tesla100GpuCount  = task_with_tesla_p100.gpuCount
    String tesla100GpuType = task_with_tesla_p100.gpuType
  }
}
