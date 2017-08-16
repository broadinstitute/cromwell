cwlVersion: "v1.0"
class: "CommandLineTool"
baseCommand: 
  - "tar"
  - "xf"
inputs: 
  - type: "string"
    inputBinding: 
      position: 2
    id: "file:///home/dan/wdl4s/tar-param.cwl#extractfile"
  - type: "File"
    inputBinding: 
      position: 1
    id: "file:///home/dan/wdl4s/tar-param.cwl#tarfile"
outputs: 
  - type: "File"
    outputBinding: 
      glob: "$(inputs.extractfile)"
    id: "file:///home/dan/wdl4s/tar-param.cwl#example_out"
id: "file:///home/dan/wdl4s/tar-param.cwl"
name: "file:///home/dan/wdl4s/tar-param.cwl"
