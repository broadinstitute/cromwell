cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: "ubuntu:latest"
baseCommand: echo
inputs:
  message:
    type: string[]?
outputs: []
