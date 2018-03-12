cwlVersion: v1.0
$graph:
- id: iwdr_array_listing
  class: CommandLineTool
  baseCommand: ['ls']
  requirements:
      - class: DockerRequirement
        dockerPull: "python:3.5.0"
      - class: InitialWorkDirRequirement
        listing: $(inputs.files)

  stdout: file_list
  inputs:
      - id: files
        type:
          type: array
          items: File
        default: [ "/Users/chrisl/Downloads/cwl/allRequirements.txt" ]
  arguments: [ "*.txt" ]
  outputs:
      - id: file_list
        type: string
        outputBinding:
          glob: file_list
          loadContents: true
          outputEval: $(self[0].contents.trim())
