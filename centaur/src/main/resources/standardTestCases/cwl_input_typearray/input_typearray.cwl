cwlVersion: v1.0
$graph:
- id: input_typearray
  cwlVersion: v1.0
  class: CommandLineTool
  baseCommand: ['/bin/echo']
  stdout: "response.txt"
  requirements:
  - class: DockerRequirement
    dockerPull: "ubuntu:latest"
  - class: InlineJavascriptRequirement
  inputs:
    value:
      type:
        - string
        - File
      inputBinding:
        position: 1
  outputs:
    response:
      type: string
      outputBinding:
        glob: response.txt
        loadContents: true
        outputEval: $(self[0].contents.split("/").slice(-1)[0])
