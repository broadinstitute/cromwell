class: Workflow
cwlVersion: v1.0
hints: []
inputs:
- id: alignment_rec
  type:
    fields:
    - name: resources
      type: string
    - name: description
      type: string
    - name: config__algorithm__align_split_size
      type:
      - 'null'
      - string
    - name: files
      type:
        items: File
        type: array
    - name: config__algorithm__trim_reads
      type:
      - string
      - 'null'
      - boolean
    - name: reference__fasta__base
      type: File
    - name: config__algorithm__adapters
      type:
      - 'null'
      - items: 'null'
        type: array
    - name: rgnames__lb
      type:
      - 'null'
      - string
    - name: rgnames__rg
      type: string
    - name: rgnames__lane
      type: string
    - name: reference__bwa__indexes
      type: File
    - name: config__algorithm__bam_clean
      type:
      - string
      - 'null'
      - boolean
    - name: config__algorithm__aligner
      type: string
    - name: rgnames__pl
      type: string
    - name: rgnames__pu
      type: string
    - name: config__algorithm__mark_duplicates
      type:
      - string
      - 'null'
      - boolean
    - name: analysis
      type: string
    - name: rgnames__sample
      type: string
    - name: config__algorithm__variant_regions
      type: File
    name: alignment_rec
    type: record
outputs:
- id: align_bam
  outputSource: merge_split_alignments/align_bam
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
- id: work_bam_plus__disc
  outputSource: merge_split_alignments/work_bam_plus__disc
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
- id: work_bam_plus__sr
  outputSource: merge_split_alignments/work_bam_plus__sr
  secondaryFiles:
  - .bai
  type:
  - File
  - 'null'
- id: hla__fastq
  outputSource: merge_split_alignments/hla__fastq
  type:
  - 'null'
  - items: File
    type: array
requirements:
- class: EnvVarRequirement
  envDef:
  - envName: MPLCONFIGDIR
    envValue: .
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
steps:
- id: prep_align_inputs
  in:
  - id: alignment_rec
    source: alignment_rec
  out:
  - id: process_alignment_rec
  run: steps/prep_align_inputs.cwl
- id: process_alignment
  in:
  - id: alignment_rec
    source: alignment_rec
  - id: process_alignment_rec
    source: prep_align_inputs/process_alignment_rec
  out:
  - id: work_bam
  - id: align_bam
  - id: hla__fastq
  - id: work_bam_plus__disc
  - id: work_bam_plus__sr
  run: steps/process_alignment.cwl
  scatter:
  - process_alignment_rec
  scatterMethod: dotproduct
- id: merge_split_alignments
  in:
  - id: alignment_rec
    source: alignment_rec
  - id: work_bam
    source: process_alignment/work_bam
  - id: align_bam_toolinput
    source: process_alignment/align_bam
  - id: work_bam_plus__disc_toolinput
    source: process_alignment/work_bam_plus__disc
  - id: work_bam_plus__sr_toolinput
    source: process_alignment/work_bam_plus__sr
  - id: hla__fastq_toolinput
    source: process_alignment/hla__fastq
  out:
  - id: align_bam
  - id: work_bam_plus__disc
  - id: work_bam_plus__sr
  - id: hla__fastq
  run: steps/merge_split_alignments.cwl
