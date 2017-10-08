cwlVersion: v1.0
class: Workflow
inputs:
- id: file:///Users/danb/wdl4s/r.cwl#pattern
  type: string
outputs:
- id: file:///Users/danb/wdl4s/r.cwl#cgrep-count
  outputSource: file:///Users/danb/wdl4s/r.cwl#cgrep/cgrep-count
  type: int
- id: file:///Users/danb/wdl4s/r.cwl#wc-count
  outputSource: file:///Users/danb/wdl4s/r.cwl#wc/wc-count
  type: int
steps:
- id: file:///Users/danb/wdl4s/r.cwl#ps
  in: []
  out:
  - file:///Users/danb/wdl4s/r.cwl#ps/ps-stdOut
  run:
    inputs: []
    outputs:
    - id: file:///Users/danb/wdl4s/r.cwl#ps/0b4ba500-5584-4fed-a831-9fa6f914ad3f/ps-stdOut
      outputBinding:
        glob: ps-stdOut.txt
      type: File
    class: CommandLineTool
    baseCommand: ps
    stdout: ps-stdOut.txt
    id: file:///Users/danb/wdl4s/r.cwl#ps/0b4ba500-5584-4fed-a831-9fa6f914ad3f
- id: file:///Users/danb/wdl4s/r.cwl#cgrep
  in:
  - id: file:///Users/danb/wdl4s/r.cwl#cgrep/pattern
    source: file:///Users/danb/wdl4s/r.cwl#pattern
  - id: file:///Users/danb/wdl4s/r.cwl#cgrep/file
    source: file:///Users/danb/wdl4s/r.cwl#ps/ps-stdOut
  out:
  - id: file:///Users/danb/wdl4s/r.cwl#cgrep/cgrep-count
  run:
    inputs:
    - id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875/pattern
      type: string
    - id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875/file
      type: File
    outputs:
    - id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875/cgrep-count
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
    id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875
- id: file:///Users/danb/wdl4s/r.cwl#wc
  in:
  - id: file:///Users/danb/wdl4s/r.cwl#wc/file
    source: file:///Users/danb/wdl4s/r.cwl#ps/ps-stdOut
  out:
  - id: file:///Users/danb/wdl4s/r.cwl#wc/wc-count
  run:
    inputs:
    - id: file:///Users/danb/wdl4s/r.cwl#wc/45d98851-7bfe-473e-ab24-aac922553f3e/file
      type: File
    outputs:
    - id: file:///Users/danb/wdl4s/r.cwl#wc/45d98851-7bfe-473e-ab24-aac922553f3e/wc-count
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
    id: file:///Users/danb/wdl4s/r.cwl#wc/45d98851-7bfe-473e-ab24-aac922553f3e
id: file:///Users/danb/wdl4s/r.cwl
