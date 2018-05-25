cwlVersion: v1.0
$graph:
- id: cwl_input_binding_expression
  class: CommandLineTool
  requirements:
    - class: InlineJavascriptRequirement
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"
  inputs:
    - id: hello
      type: string
      inputBinding:
        valueFrom: $(self + " world")
        position: 1
  outputs:
    b:
      type: string
      outputBinding:
        loadContents: true
        glob: stdout
        outputEval: $(self[0].contents.trim())
  stdout: stdout
  baseCommand: []
  arguments: [echo]
