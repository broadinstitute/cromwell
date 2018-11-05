class: Workflow
cwlVersion: v1.0
hints: []
inputs:
- id: files
  type:
    items:
      items: File
      type: array
    type: array
- id: reference__fasta__base
  secondaryFiles:
  - .fai
  - ^.dict
  type:
    items: File
    type: array
- id: config__algorithm__expression_caller
  type:
    items:
      items: string
      type: array
    type: array
- id: resources
  type:
    items: string
    type: array
- id: config__algorithm__coverage_interval
  type:
    items:
    - 'null'
    - string
    type: array
- id: rgnames__lb
  type:
    items:
    - 'null'
    - string
    type: array
- id: rgnames__rg
  type:
    items: string
    type: array
- id: reference__hisat2__indexes
  type:
    items: File
    type: array
- id: config__algorithm__fusion_caller
  type:
    items:
      items: string
      type: array
    type: array
- id: description
  type:
    items: string
    type: array
- id: config__algorithm__aligner
  type:
    items: string
    type: array
- id: rgnames__pl
  type:
    items: string
    type: array
- id: genome_build
  type:
    items: string
    type: array
- id: rgnames__pu
  type:
    items: string
    type: array
- id: genome_resources__rnaseq__transcripts
  secondaryFiles:
  - .db
  type:
    items: File
    type: array
- id: config__algorithm__quality_format
  type:
    items: string
    type: array
- id: config__algorithm__tools_off
  type:
    items:
    - 'null'
    - items: 'null'
      type: array
    type: array
- id: config__algorithm__qc
  type:
    items:
      items: string
      type: array
    type: array
- id: analysis
  type:
    items: string
    type: array
- id: rgnames__sample
  type:
    items: string
    type: array
- id: config__algorithm__tools_on
  type:
    items:
    - 'null'
    - items: 'null'
      type: array
    type: array
- id: rgnames__lane
  type:
    items: string
    type: array
outputs:
- id: rgnames__sample_out
  outputSource: rgnames__sample
  type:
    items: string
    type: array
- id: align_bam
  outputSource: process_alignment/align_bam
  type:
    items: File
    type: array
- id: quant__tsv
  outputSource: rnaseq_quantitate/quant__tsv
  type:
    items: File
    type: array
- id: summary__multiqc
  outputSource: multiqc_summary/summary__multiqc
  type:
    items:
    - File
    - 'null'
    type: array
requirements:
- class: EnvVarRequirement
  envDef:
  - envName: MPLCONFIGDIR
    envValue: .
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
steps:
- id: prepare_sample
  in:
  - id: files
    source: files
  - id: rgnames__sample
    source: rgnames__sample
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: genome_build
    source: genome_build
  - id: genome_resources__rnaseq__transcripts
    source: genome_resources__rnaseq__transcripts
  - id: analysis
    source: analysis
  - id: rgnames__pl
    source: rgnames__pl
  - id: rgnames__pu
    source: rgnames__pu
  - id: rgnames__lane
    source: rgnames__lane
  - id: rgnames__rg
    source: rgnames__rg
  - id: rgnames__lb
    source: rgnames__lb
  - id: reference__hisat2__indexes
    source: reference__hisat2__indexes
  - id: config__algorithm__aligner
    source: config__algorithm__aligner
  - id: config__algorithm__expression_caller
    source: config__algorithm__expression_caller
  - id: config__algorithm__fusion_caller
    source: config__algorithm__fusion_caller
  - id: config__algorithm__quality_format
    source: config__algorithm__quality_format
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: prep_rec
  run: steps/prepare_sample.cwl
  scatter:
  - files
  - rgnames__sample
  - reference__fasta__base
  - genome_build
  - genome_resources__rnaseq__transcripts
  - analysis
  - rgnames__pl
  - rgnames__pu
  - rgnames__lane
  - rgnames__rg
  - rgnames__lb
  - reference__hisat2__indexes
  - config__algorithm__aligner
  - config__algorithm__expression_caller
  - config__algorithm__fusion_caller
  - config__algorithm__quality_format
  - resources
  - description
  scatterMethod: dotproduct
- id: trim_sample
  in:
  - id: prep_rec
    source: prepare_sample/prep_rec
  out:
  - id: trim_rec
  run: steps/trim_sample.cwl
  scatter:
  - prep_rec
  scatterMethod: dotproduct
- id: process_alignment
  in:
  - id: trim_rec
    source: trim_sample/trim_rec
  out:
  - id: align_bam
  run: steps/process_alignment.cwl
  scatter:
  - trim_rec
  scatterMethod: dotproduct
- id: rnaseq_quantitate
  in:
  - id: trim_rec
    source: trim_sample/trim_rec
  - id: align_bam
    source: process_alignment/align_bam
  out:
  - id: count_file
  - id: quant__tsv
  - id: quant__hdf5
  - id: quant__fusion
  run: steps/rnaseq_quantitate.cwl
  scatter:
  - trim_rec
  - align_bam
  scatterMethod: dotproduct
- id: qc_to_rec
  in:
  - id: align_bam
    source: process_alignment/align_bam
  - id: analysis
    source: analysis
  - id: reference__fasta__base
    source: reference__fasta__base
  - id: genome_resources__rnaseq__transcripts
    source: genome_resources__rnaseq__transcripts
  - id: genome_build
    source: genome_build
  - id: config__algorithm__coverage_interval
    source: config__algorithm__coverage_interval
  - id: config__algorithm__tools_on
    source: config__algorithm__tools_on
  - id: config__algorithm__tools_off
    source: config__algorithm__tools_off
  - id: config__algorithm__qc
    source: config__algorithm__qc
  - id: resources
    source: resources
  - id: description
    source: description
  out:
  - id: qc_rec
  run: steps/qc_to_rec.cwl
- id: pipeline_summary
  in:
  - id: qc_rec
    source: qc_to_rec/qc_rec
  out:
  - id: qcout_rec
  run: steps/pipeline_summary.cwl
  scatter:
  - qc_rec
  scatterMethod: dotproduct
- id: multiqc_summary
  in:
  - id: qcout_rec
    source: pipeline_summary/qcout_rec
  out:
  - id: summary__multiqc
  run: steps/multiqc_summary.cwl
- id: detect_fusions
  in:
  - id: quant__fusion
    source: rnaseq_quantitate/quant__fusion
  - id: quant__hdf5
    source: rnaseq_quantitate/quant__hdf5
  - id: trim_rec
    source: trim_sample/trim_rec
  out:
  - id: fusion__fasta
  - id: fusion__json
  run: steps/detect_fusions.cwl
  scatter:
  - quant__fusion
  - quant__hdf5
  - trim_rec
  scatterMethod: dotproduct
