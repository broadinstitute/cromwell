#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

inputs:
- id: pin
  type: string
  default: "v"

outputs:
    count_output:
      type: int
      outputSource: threestep/wc-count

requirements:
  - class: SubworkflowFeatureRequirement

steps:
  threestep:
    run: three_step.cwl
    in:
    - id: pattern
      source: pin
    out: [wc-count]
