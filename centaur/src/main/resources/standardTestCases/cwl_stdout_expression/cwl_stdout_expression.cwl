cwlVersion: v1.0
$graph:
- id: cwl_stdout_expression
  class: CommandLineTool
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"
  inputs:
    - id: foo
      type: string
    - id: bar
      type: string
  outputs:
    b: stdout
  stdout: "stdout-$(inputs.foo)-$(inputs.bar).txt"
  baseCommand: []
  arguments: [echo, $(inputs.foo), $(inputs.bar)]
