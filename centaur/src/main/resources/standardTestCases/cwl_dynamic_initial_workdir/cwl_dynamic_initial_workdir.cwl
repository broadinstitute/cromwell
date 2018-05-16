cwlVersion: v1.0
$graph:
  id: main
  class: CommandLineTool
  cwlVersion: v1.0
  requirements:
    - class: ShellCommandRequirement
    - class: InitialWorkDirRequirement
      listing: $(inputs.indir.listing)
    - class: DockerRequirement
      dockerPull: "ubuntu:latest"
  inputs:
    indir: Directory
  outputs:
    sha:
      type: string
      outputBinding:
        glob: output.txt
        loadContents: true
        outputEval: $(self[0].contents)
  arguments: ["find", "-L", ".", "-name", "?",
    {shellQuote: false, valueFrom: "|"},
    "sort",
    {shellQuote: false, valueFrom: "|"},
    "sha1sum",
    {shellQuote: false, valueFrom: "|"},
    "cut", "-c", "1-40"
  ]
  stdout: output.txt
