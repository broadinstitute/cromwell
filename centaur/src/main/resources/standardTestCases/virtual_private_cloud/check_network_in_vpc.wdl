version 1.0

task get_network {
  command {
    set -euo pipefail

    apt-get install --assume-yes jq > /dev/null
    INSTANCE=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
    ZONE=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/zone" -H "Metadata-Flavor: Google" | sed -E 's!.*/(.*)!\1!')
    TOKEN=$(gcloud auth application-default print-access-token)
    INSTANCE_METADATA=$(curl "https://www.googleapis.com/compute/v1/projects/broad-dsde-cromwell-dev/zones/$ZONE/instances/$INSTANCE" -H "Authorization: Bearer $TOKEN" -H 'Accept: application/json')
    NETWORK_OBJECT=$(echo $INSTANCE_METADATA | jq --raw-output --exit-status '.networkInterfaces[0]')
    echo $NETWORK_OBJECT | jq --exit-status '.network' | sed -E 's!.*/(.*)!\1!' | sed 's/"//g' > network
    echo $NETWORK_OBJECT | jq --exit-status '.subnetwork' | sed -E 's!.*/(.*)!\1!' | sed 's/"//g' > subnetwork
    echo $ZONE > zone
  }

  runtime {
    docker: "google/cloud-sdk:slim"
    backend: "Papiv2-Virtual-Private-Cloud"
  }

  output {
    String networkName = read_string("network")
    String subnetworkName = read_string("subnetwork")
    String zone = read_string("zone")
  }
}

workflow check_network_in_vpc {
  call get_network

   output {
     String network_used = get_network.networkName
     String subnetwork_used = get_network.subnetworkName
     String zone_used = get_network.zone
   }
}

