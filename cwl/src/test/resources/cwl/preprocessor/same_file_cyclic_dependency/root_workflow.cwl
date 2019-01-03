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
      run: "#echo-workflow-1"
      in:
        tool: tool
        in:
          valueFrom: $(inputs.tool.nameroot)
      out: [out]

- id: echo-workflow-1
  class: Workflow
  requirements:
    - class: StepInputExpressionRequirement
  inputs:
    tool: File
  outputs: []
  steps:
    root:
      run: "#echo-workflow-2"
      in:
        tool: tool
        in:
          valueFrom: $(inputs.tool.nameroot)
      out: [out]
      