cwlVersion: v1.0
class: Workflow
requirements:
  - class: StepInputExpressionRequirement
inputs:
  tool: File
outputs: []
steps:
  root:
    run: root_workflow.cwl
    in:
      tool: tool
      in:
        valueFrom: $(inputs.tool.nameroot)
    out: [out]
    