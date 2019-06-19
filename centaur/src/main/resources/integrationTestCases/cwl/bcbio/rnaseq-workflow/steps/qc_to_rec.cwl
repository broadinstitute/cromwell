$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=qc_rec:resources;description;reference__fasta__base;config__algorithm__coverage_interval;genome_build;genome_resources__rnaseq__transcripts;config__algorithm__tools_off;config__algorithm__qc;analysis;config__algorithm__tools_on;align_bam
- sentinel_inputs=align_bam:var,analysis:var,reference__fasta__base:var,genome_resources__rnaseq__transcripts:var,genome_build:var,config__algorithm__coverage_interval:var,config__algorithm__tools_on:var,config__algorithm__tools_off:var,config__algorithm__qc:var,resources:var,description:var
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- qc_to_rec
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 1
  outdirMin: 1028
  ramMin: 2048
  tmpdirMin: 2
- class: dx:InputResourceRequirement
  indirMin: 0
inputs:
- id: align_bam
  secondaryFiles:
  - .bai
  type:
    items: File
    type: array
- id: analysis
  type:
    items: string
    type: array
- id: reference__fasta__base
  secondaryFiles:
  - .fai
  - ^.dict
  type:
    items: File
    type: array
- id: genome_resources__rnaseq__transcripts
  secondaryFiles:
  - .db
  type:
    items: File
    type: array
- id: genome_build
  type:
    items: string
    type: array
- id: config__algorithm__coverage_interval
  type:
    items:
    - 'null'
    - string
    type: array
- id: config__algorithm__tools_on
  type:
    items:
    - 'null'
    - items: 'null'
      type: array
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
- id: resources
  type:
    items: string
    type: array
- id: description
  type:
    items: string
    type: array
outputs:
- id: qc_rec
  type:
    items:
      fields:
      - name: resources
        type: string
      - name: description
        type: string
      - name: reference__fasta__base
        type: File
      - name: config__algorithm__coverage_interval
        type:
        - 'null'
        - string
      - name: genome_build
        type: string
      - name: genome_resources__rnaseq__transcripts
        type: File
      - name: config__algorithm__tools_off
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: config__algorithm__qc
        type:
          items: string
          type: array
      - name: analysis
        type: string
      - name: config__algorithm__tools_on
        type:
        - 'null'
        - items: 'null'
          type: array
      - name: align_bam
        type: File
      name: qc_rec
      type: record
    type: array
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
