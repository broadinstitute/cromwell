cwlVersion: v1.0
$graph:  
- id: outputJson
  # slightly modified version of conformance test #2 to remove the default File that doesn't work on PAPI
  # It is instead being passed as an input through the input json file
  class: CommandLineTool
  cwlVersion: v1.0
  hints:
    - class: DockerRequirement
      dockerPull: python:2-slim
  inputs:
    - id: reference
      type: File
      inputBinding: { position: 2 }
  
    - id: reads
      type:
        type: array
        items: File
        inputBinding: { prefix: "-YYY" }
      inputBinding: { position: 3, prefix: "-XXX" }
  
    - id: "args.py"
      type: File
      inputBinding:
        position: -1
  
  outputs:
    # note the absence of any sort of valueFrom.
    # The output value is generated from the "cwl.output.json" file created by the python script
    args: string[]
  
  baseCommand: python
  arguments: ["bwa", "mem"]
