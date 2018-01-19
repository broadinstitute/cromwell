cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: "ubuntu:latest"
baseCommand: echo
inputs:
  message:
    type: string[]?
  unsupplied_optional:
    type: string[]?
outputs: []
