version 1.0

task get_network_labels_backend {
  meta { volatile: true }
  input {
    Array[String] commandScript
    String dockerImage
  }
  command <<<~{sep="\n" commandScript}>>>
  runtime {
    docker: dockerImage
    backend: "Papiv2-Virtual-Private-Cloud-Labels"
  }
  output {
    String networkName = read_string("network")
    String subnetworkName = read_string("subnetwork")
    String zone = read_string("zone")
  }
}

task get_network_literals_backend {
  meta { volatile: true }
  input {
    Array[String] commandScript
    String dockerImage
  }
  command <<<~{sep="\n" commandScript}>>>
  runtime {
    docker: dockerImage
    backend: "Papiv2-Virtual-Private-Cloud-Literals"
  }
  output {
    String networkName = read_string("network")
    String subnetworkName = read_string("subnetwork")
    String zone = read_string("zone")
  }
}

workflow check_network_in_vpc {
  # Create a reusable script for multiple backends using workarounds for Cromwell 66:
  #  - can't pass `backend` as a variable runtime attribute; the call will run on the default backend
  #  - use an `Array[String]` to simulate a WDL multiline string
  #  - an escaped backslash `\\` at the end of a WDL string returns a syntax error (maybe from womtool describe?)
  String backslash = "\u005c"
  Array[String] commandScript = [
    "set -euo pipefail",
    "",
    "apt-get install --assume-yes jq > /dev/null",
    "PROJECT=$(",
    "  curl " + backslash,
    "    -s \"http://metadata.google.internal/computeMetadata/v1/project/project-id\" " + backslash,
    "    -H \"Metadata-Flavor: Google\" |",
    "  sed -E 's!.*/(.*)!\\1!'",
    ")",
    "ZONE=$(",
    "  curl " + backslash,
    "    -s \"http://metadata.google.internal/computeMetadata/v1/instance/zone\" " + backslash,
    "    -H \"Metadata-Flavor: Google\" |",
    "  sed -E 's!.*/(.*)!\\1!'",
    ")",
    "INSTANCE=$(",
    "  curl " + backslash,
    "    -s \"http://metadata.google.internal/computeMetadata/v1/instance/name\" " + backslash,
    "    -H \"Metadata-Flavor: Google\"",
    ")",
    "TOKEN=$(gcloud auth application-default print-access-token)",
    "INSTANCE_METADATA=$(",
    "  curl \"https://www.googleapis.com/compute/v1/projects/$PROJECT/zones/$ZONE/instances/$INSTANCE\" " + backslash,
    "    -H \"Authorization: Bearer $TOKEN\" " + backslash,
    "    -H 'Accept: application/json'",
    ")",
    "NETWORK_OBJECT=$(echo $INSTANCE_METADATA | jq --raw-output --exit-status '.networkInterfaces[0]')",
    "echo $NETWORK_OBJECT | jq --exit-status '.network' | sed -E 's!.*/(.*)!\\1!' | sed 's/\"//g' > network",
    "echo $NETWORK_OBJECT | jq --exit-status '.subnetwork' | sed -E 's!.*/(.*)!\\1!' | sed 's/\"//g' > subnetwork",
    "echo $ZONE > zone"
  ]
  String dockerImage = "gcr.io/google.com/cloudsdktool/cloud-sdk:slim"

  call get_network_labels_backend {
    input:
      commandScript = commandScript,
      dockerImage = dockerImage
  }

  call get_network_literals_backend {
    input:
      commandScript = commandScript,
      dockerImage = dockerImage
  }

   output {
     String network_used_labels = get_network_labels_backend.networkName
     String subnetwork_used_labels = get_network_labels_backend.subnetworkName
     String zone_used_labels = get_network_labels_backend.zone
     String network_used_literals = get_network_literals_backend.networkName
     String subnetwork_used_literals = get_network_literals_backend.subnetworkName
     String zone_used_literals = get_network_literals_backend.zone
   }
}

