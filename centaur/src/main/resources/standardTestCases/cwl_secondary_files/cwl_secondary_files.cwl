cwlVersion: v1.0
$graph:
- id: main
  class: CommandLineTool
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"
  inputs:
    - id: command
      type: string
    - id: f
      type: File
      inputBinding:
        position: 2
      secondaryFiles: [.also]
    - id: fs
      type:
        type: array
        items: File
        inputBinding:
          position: 3
      secondaryFiles: [.also]
  outputs:
    the_answer:
      type: string
      outputBinding:
        outputEval: ${ return "$(" + 42 + ")"; }
  baseCommand: []
  arguments: ["bash", "-c", $(inputs.command)]
