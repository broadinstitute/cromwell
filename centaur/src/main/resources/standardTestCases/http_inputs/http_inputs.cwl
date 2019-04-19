cwlVersion: v1.0
$graph:
  - id: http_inputs
    class: Workflow
    inputs:
    - id: jamie
      type: File
    outputs:
    - id: md5
      outputSource: "#http_inputs/sum/md5"
      type: int
    steps:
    - id: sum
      in:
      - id: jamie
        source: "#http_inputs/jamie"
      out:
      - id: md5
      # Workflow step-level DockerRequirement
      requirements:
        DockerRequirement:
          dockerPull: "ubuntu:latest"
      run:
        inputs:
        - id: jamie
          type: File
        outputs:
        - id: md5
          outputBinding:
            glob: md5-stdOut.txt
            loadContents: true
            outputEval: $(self[0].contents)
          type: string
        class: CommandLineTool
        requirements:
        - class: ShellCommandRequirement
        - class: InlineJavascriptRequirement
        arguments:
        - valueFrom: "/usr/bin/md5sum"
          shellQuote: false
        - valueFrom: $(inputs.jamie)
          shellQuote: true
        - valueFrom: '|'
          shellQuote: false
        - valueFrom: cut
          shellQuote: false
        - valueFrom: -c1-32
          shellQuote: false
        stdout: md5-stdOut.txt
