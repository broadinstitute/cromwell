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
  arguments:
      - position: 3
        valueFrom: "sentinel"
  inputs:
    value_f:
      type:
        - string
        - File
      inputBinding:
        position: 1
      doc: "an input to test with a File value"
    value_s:
      type:
        - string
        - File
      inputBinding:
        position: 2
      doc: "an input to test with a string value"
  outputs:
    response_f:
      type: string
      outputBinding:
        glob: response.txt
        loadContents: true
        outputEval: $(self[0].contents.split(" ")[0].split("/").slice(-1)[0])
    response_s:
      type: string
      outputBinding:
        glob: response.txt
        loadContents: true
        outputEval: $(self[0].contents.split(" ")[1].split("/").slice(-1)[0])
