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
# size is optional according to the spec, however there are real world CWLs using it
arguments: ["echo", "$(inputs.infile.size)"]
stdout: output.txt
