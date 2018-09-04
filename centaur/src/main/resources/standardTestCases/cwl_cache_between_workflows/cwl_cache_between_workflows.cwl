cwlVersion: v1.0
$graph:
- id: average
  class: CommandLineTool

  inputs:
    base1:
      type: int
      default: 9
    base2:
      type: int
      default: 13

  outputs:
    baseAverage:
      type: float
      outputBinding:
        glob: stdout.txt
        loadContents: true
        outputEval: $(parseFloat(self[0].contents))

  arguments:
  - valueFrom: "echo $((inputs.base1 + inputs.base2) / 2)"
    shellQuote: false

  stdout: stdout.txt

- id: height-product
  class: CommandLineTool

  inputs:
    baseAverage: float
    height:
      type: float
      default: 7.0

  outputs:
    trapezoidalArea:
      type: long
      outputBinding:
        glob: stdout.txt
        loadContents: true
        outputEval: $(parseFloat(self[0].contents))

  arguments:
  - valueFrom: "echo $(inputs.baseAverage * inputs.height)"
    shellQuote: false

  stdout: stdout.txt

- id: cwl-cache-between-workflows
  class: Workflow

  requirements:
  - class: DockerRequirement
    dockerPull: "ubuntu:latest"
  - class: ShellCommandRequirement

  inputs: []

  outputs:
    baseAverage:
      type: float
      outputSource: step-average/baseAverage
    trapezoidalArea:
      type: long
      outputSource: step-product/trapezoidalArea

  steps:
    step-average:
      run: "#average"
      in: []
      out: [baseAverage]

    step-product:
      run: "#height-product"
      in:
        baseAverage: step-average/baseAverage
      out: [trapezoidalArea]
