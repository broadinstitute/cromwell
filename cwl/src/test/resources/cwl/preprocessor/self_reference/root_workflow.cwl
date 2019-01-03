cwlVersion: v1.0
$graph:
- id: echo-workflow-2
  class: Workflow
  requirements:
    - class: StepInputExpressionRequirement
  inputs:
    tool: File
  outputs: []
  steps:
    root:
      run: "#echo-tool"
      in:
        tool: tool
        in:
          valueFrom: $(inputs.tool.nameroot)
      out: [out]

- id: echo-tool
  class: CommandLineTool
  baseCommand: [echo]
  inputs:
    in:
      type: string
      inputBinding:
        position: 1
  outputs:
    out:
      type: string
      valueFrom: "hello"
      