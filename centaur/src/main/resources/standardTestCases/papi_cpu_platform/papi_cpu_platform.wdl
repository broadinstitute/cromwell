version 1.0

task cpu_platform {
    input {
        String cpu_platform
    }
    command {
       apt-get install --assume-yes jq > /dev/null
       NAME=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/name)
       FULLY_QUALIFIED_ZONE=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone)
       ZONE=$(basename "$FULLY_QUALIFIED_ZONE")
       PROJECT=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/project/project-id)
       curl -s -H "Authorization: Bearer $(gcloud auth print-access-token)" "https://www.googleapis.com/compute/v1/projects/$PROJECT/zones/$ZONE/instances/$NAME?fields=cpuPlatform" | jq -r '.cpuPlatform'
    }

    runtime {
        docker: "google/cloud-sdk:slim"
        cpuPlatform: cpu_platform
    }

    output {
        String cpuPlatform = read_string(stdout())
    }
}

workflow cpus {
   call cpu_platform as haswell { input: cpu_platform = "Intel Haswell" }
   call cpu_platform as broadwell { input: cpu_platform = "Intel Broadwell" }
}
