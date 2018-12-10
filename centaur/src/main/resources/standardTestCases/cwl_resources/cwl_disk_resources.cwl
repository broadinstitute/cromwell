cwlVersion: v1.0
# Asks for 10 GB (~=10240 MB) of output dir, tmp dir and input dir.
# Unfortunately the internal metadata from the VM, despite having some information about the disks doesn't give their size.
# To avoid possibly unreliable computations based on free or other tool, instead use the GCE API directly to look up the disk info.
# To do that we need project name, zone and name of the disk which we all obtain via the internal metadata
# Of importance is the assumption that <disk name> == <instance name>-1
# Note the "-1" after the instance name. That is because the disk just named <instance name> is the boot disk. The working disk has a "-1" suffix
# The test expectation ensures that the output is 30, which means all 3 resources have been taken into account.
$namespaces:
  dx: https://www.dnanexus.com/cwl#
$graph:  
- id: diskSizeTool
  class: CommandLineTool
  cwlVersion: v1.0
  doc: "Asks for disk minimums"
  requirements:
    ResourceRequirement:
      outdirMin: 10240
      tmpdirMin: 10240
  hints:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: google/cloud-sdk:slim
  - class: dx:InputResourceRequirement
    indirMin: 10240
  inputs: []
  outputs:
    disk_size:
      type: int
      outputBinding:
        glob: stdout
        loadContents: true
        outputEval: $(parseInt(self[0].contents.trim()))
  arguments:
    - valueFrom: >
        apt-get install --assume-yes jq > /dev/null && NAME=`curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/name` && ZONE=`basename \`curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone\`` && PROJECT=`curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/project/project-id` && curl -s -H "Authorization: Bearer `gcloud auth print-access-token`" "https://www.googleapis.com/compute/v1/projects/${PROJECT}/zones/${ZONE}/disks/${NAME}-1?fields=sizeGb" | jq -r '.sizeGb'
      shellQuote: false
  stdout: stdout
