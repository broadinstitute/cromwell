#!/usr/bin/env cwl-runner

cwlVersion: v1.0

requirements:
  - class: DockerRequirement
    dockerPull: ubuntu:bionic-20180426

class: CommandLineTool

inputs:
  - id: input
    type: string
    inputBinding:
      position: 0

outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.input)
      
baseCommand: [touch]
