cwlVersion: v1.0
class: Workflow
id: three_step
# Workflow-level DockerRequirement
hints:
  DockerRequirement:
    dockerPull: "ubuntu:latest"
inputs:
- id: pattern
  type: string
outputs:
- id: cgrep-count
  outputSource: "#cgrep/cgrep-count"
  type: int
- id: wc-count
  outputSource: "#wc/wc-count"
  type: int
steps:
- id: ps
  in: []
  out:
  - id: ps-stdOut
  run:
    inputs: []
    outputs:
    - id: ps-stdOut
      outputBinding:
        glob: ps-stdOut.txt
      type: File
    class: CommandLineTool
    requirements:
      # Command line tool-level DockerRequirement
      # Check it: this DockerRequirement does not use the `class` formatting but should still parse
      # thanks to the magic of SALAD.
      DockerRequirement:
        dockerPull: "ubuntu:bionic"
    baseCommand: ps
    stdout: ps-stdOut.txt
- id: cgrep
  in:
  - id: pattern
    source: "#pattern"
  - id: file
    source: "#ps/ps-stdOut"
  out:
  - id: cgrep-count
  # Workflow step-level DockerRequirement
  requirements:
    DockerRequirement:
      dockerPull: "debian:jessie"
  run:
    inputs:
    - id: pattern
      type: string
    - id: file
      type: File
    outputs:
    - id: cgrep-count
      outputBinding:
        glob: cgrep-stdOut.txt
        loadContents: true
        outputEval: $(parseInt(self[0].contents))
      type: int
    class: CommandLineTool
    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    arguments:
    - valueFrom: grep
      shellQuote: false
    - valueFrom: $(inputs.pattern)
      shellQuote: false
    - valueFrom: ${return inputs.file}
      shellQuote: false
    - valueFrom: '|'
      shellQuote: false
    - valueFrom: wc
      shellQuote: false
    - valueFrom: -l
      shellQuote: false
    stdout: cgrep-stdOut.txt
- id: wc
  in:
  - id: file
    source: "#ps/ps-stdOut"
  out:
  - id: wc-count
  run:
    inputs:
    - id: file
      type: File
    outputs:
    - id: wc-count
      outputBinding:
        glob: wc-stdOut.txt
        loadContents: true
        outputEval: $(parseInt(self[0].contents))
      type: int
    class: CommandLineTool
    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    arguments:
    - valueFrom: cat
      shellQuote: false
    - valueFrom: $(inputs.file)
      shellQuote: false
    - valueFrom: '|'
      shellQuote: false
    - valueFrom: wc
      shellQuote: false
    - valueFrom: -l
      shellQuote: false
    stdout: wc-stdOut.txt
