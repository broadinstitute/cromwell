cwlVersion: v1.0
$graph:
- id: interpolatedStrings
  class: CommandLineTool
  cwlVersion: v1.0
  requirements:
    - class: ShellCommandRequirement
  hints:
    DockerRequirement:
      dockerPull: "debian:stretch-slim"

  inputs:
    - id: bar
      type: string
    - id: qux
      type: string

  outputs:
    - id: rfc3092
      type: string

  arguments:
     - valueFrom: >
         echo '{ "rfc3092": "foo $(inputs.bar) baz $(inputs.qux) quux" }' > cwl.output.json
       shellQuote: false
