class: Workflow
cwlVersion: v1.0
inputs:
  file1: File
outputs:
  count_output: {type: int, outputSource: step0/count_output}
requirements:
  SubworkflowFeatureRequirement: {}
steps:
  step0:
    in: {file1: file1}
    out: [count_output]
    run:
      class: Workflow
      inputs:
        file1: File
      outputs:
        count_output: {type: int, outputSource: step2/output}
      steps:
        step1: {run: wc-tool.cwl, in: {file1: file1}, out: [output]}
        step2: {run: parseInt-tool.cwl, in: {file1: step1/output}, out: [output]}
