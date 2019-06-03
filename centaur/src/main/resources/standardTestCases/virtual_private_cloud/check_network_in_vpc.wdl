version 1.0

task get_network {
  command {
    apt-get install --assume-yes jq > /dev/null
    INSTANCE=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
    ZONE=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/zone" -H "Metadata-Flavor: Google" | sed -E 's!.*/(.*)!\1!')
    TOKEN=$(gcloud auth application-default print-access-token)
    INSTANCE_METADATA=$(curl "https://www.googleapis.com/compute/v1/projects/broad-dsde-cromwell-dev/zones/$ZONE/instances/$INSTANCE" -H "Authorization: Bearer $TOKEN" -H 'Accept: application/json')
    echo $INSTANCE_METADATA | jq -r '.networkInterfaces[0].network' | sed -E 's!.*/(.*)!\1!'
  }

  runtime {
    docker: "google/cloud-sdk:slim"
    backend: "Papiv2-Virtual-Private-Cloud"
  }

  output {
    String network = read_string(stdout())
  }
}

workflow check_network_in_vpc {
  call get_network

   output {
     String network_used = get_network.network
   }
}
