class: CommandLineTool
cwlVersion: v1.0
requirements:
  - class: ShellCommandRequirement
inputs:
  - id: infile
    type: File
    default:
      class: File
      location: exampledir/example.txt
outputs:
  outlist:
    type: File
    outputBinding:
      glob: output.txt
arguments: ["echo", "$(inputs.infile.basename)"]
stdout: output.txt
