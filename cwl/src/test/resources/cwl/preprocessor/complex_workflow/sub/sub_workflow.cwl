cwlVersion: v1.0
$graph:
- id: sub-echo-workflow-2
  class: Workflow
  requirements:
    - class: StepInputExpressionRequirement
  inputs:
    tool: File
  outputs: []
  steps:
    root:
      run: ../../echo_tool.cwl
      in:
        tool: tool
        in:
          valueFrom: $(inputs.tool.nameroot)
      out: [out]

- id: sub-echo-workflow-1
  class: Workflow
  requirements:
    - class: StepInputExpressionRequirement
  inputs:
    tool: File
  outputs: []
  steps:
    root:
      run: "../root_workflow.cwl#echo-workflow-1"
      in:
        tool: tool
        in:
          valueFrom: $(inputs.tool.nameroot)
      out: [out]
      