cwlVersion: "v1.0"
class: "CommandLineTool"
baseCommand: ["tar", "xf"]
inputs: 
  - type: "string"
    inputBinding: 
      position: 2
    id: "extractfile"
  - type: "File"
    inputBinding: 
      position: 1
    id: "tarfile"
outputs: 
  - type: "File"
    outputBinding: 
      glob: "$(inputs.extractfile)"
    id: "example_out"
