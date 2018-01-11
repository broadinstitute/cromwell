class: CommandLineTool
cwlVersion: v1.0
requirements:
  - class: ShellCommandRequirement
inputs:
  - id: indir
    type: Directory
    default:
      class: Directory
      location: exampledir
outputs:
  outlist:
    type: File
    outputBinding:
      glob: output.txt
arguments: ["echo", "$(inputs.indir.basename)"]
stdout: output.txt
