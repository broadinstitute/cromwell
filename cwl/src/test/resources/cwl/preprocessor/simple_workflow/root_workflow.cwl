cwlVersion: v1.0
class: Workflow
requirements:
  - class: StepInputExpressionRequirement
inputs:
   tool: File
outputs: []
steps:
  step1:
    run: ../echo_tool.cwl
    in:
      tool: tool
      in:
        valueFrom: $(inputs.tool.nameroot)
    out: [out]
  step2:
    run: ../echo_tool.cwl
    in:
      tool: tool
      in:
        valueFrom: $(inputs.tool.nameroot)
    out: [out]
    