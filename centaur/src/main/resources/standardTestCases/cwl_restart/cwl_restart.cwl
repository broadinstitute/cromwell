cwlVersion: v1.0
class: Workflow
id: cwl_restart
# Workflow-level DockerRequirement
requirements:
  DockerRequirement:
    dockerPull: "ubuntu:latest"
inputs: []
outputs:
- id: flag
  outputSource: "#third_step/flag"
  type: boolean
steps:
- id: first_step
  in: []
  out:
  - id: firstfoofile
  # A dash is not valid in a WDL identifier and would have broken older `IdentifierAndPathPattern` parsing.
  - id: first-bar-file
  run:
    inputs: []
    outputs:
    - id: firstfoofile
      outputBinding:
        glob: out.txt
      type: File
    - id: first-bar-file
      outputBinding:
        glob: out.txt
      type: File
    class: CommandLineTool
    requirements:
    - class: ShellCommandRequirement
    arguments:
    - valueFrom: echo
    - valueFrom: "I have a bad feeling about this"
    - valueFrom: ">"
      shellQuote: false
    - valueFrom: "out.txt"
- id: cromwell_killer
  in:
  - id: fooinput
    source: "#first_step/firstfoofile"
  - id: bar-input
    source: "#first_step/first-bar-file"
  out:
  - id: "footxt"
  - id: "bar-txt"
  run:
    inputs:
    - id: fooinput
      type: File
    - id: bar-input
      type: File
    outputs:
    - id: "footxt"
      outputBinding:
        glob: foo.txt
      type: File
    - id: "bar-txt"
      outputBinding:
        glob: bar.txt
      type: File
    class: CommandLineTool
    requirements:
    - class: ShellCommandRequirement
    arguments:
    - valueFrom: sleep
    - valueFrom: "60"
    - valueFrom: "&&"
      shellQuote: false
    - valueFrom: echo
    - valueFrom: foo
    - valueFrom: ">"
      shellQuote: false
    - valueFrom: "foo.txt"
    - valueFrom: "&&"
      shellQuote: false
    - valueFrom: echo
    - valueFrom: bar
    - valueFrom: ">"
      shellQuote: false
    - valueFrom: "bar.txt"
- id: third_step
  in:
  - id: foo
    source: "#cromwell_killer/footxt"
  - id: bar
    source: "#cromwell_killer/bar-txt"
  out:
  - id: flag
  run:
    inputs:
    - id: foo
      type: File
    - id: bar
      type: File
    outputs:
    - id: flag
      outputBinding:
        outputEval: $(true)
      type: boolean
    class: CommandLineTool
    requirements:
    - class: InlineJavascriptRequirement
    - class: ShellCommandRequirement
    arguments:
    - valueFrom: echo
    - valueFrom: "Are we alive???"
    - valueFrom: "&&"
      shellQuote: false
    - valueFrom: echo
    - valueFrom: $(inputs.foo)
    - valueFrom: "&&"
      shellQuote: false
    - valueFrom: echo
    - valueFrom: $(inputs.bar)
