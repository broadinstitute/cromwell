cwlVersion: v1.0
$graph:  
- id: machineTypeTool
  class: CommandLineTool
  cwlVersion: v1.0
  doc: "Asks for CPU and Memory minimums"
  requirements:
    ResourceRequirement:
      coresMin: 2
      coresMax: 2
      ramMin: 7GB
      ramMax: 7GB
  hints:
    DockerRequirement:
      dockerPull: python:latest
  inputs: []
  outputs:
    machine_type:
      type: string
      outputBinding:
        glob: stdout
        loadContents: true
        outputEval: $(self[0].contents.trim())
  baseCommand: ['curl', 'http://metadata.google.internal/computeMetadata/v1/instance/machine-type', '-H', 'Metadata-Flavor: Google']
  stdout: stdout


