$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=fusion__fasta,fusion__json
- sentinel_inputs=quant__fusion:var,quant__hdf5:var,trim_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- detect_fusions
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 2
  outdirMin: 1024
  ramMin: 4096
  tmpdirMin: 0
- class: dx:InputResourceRequirement
  indirMin: 4
- class: SoftwareRequirement
  packages:
  - package: pizzly
    specs:
    - https://anaconda.org/bioconda/pizzly
inputs:
- id: quant__fusion
  type: File
- id: quant__hdf5
  type: File
- id: trim_rec
  type:
    fields:
    - name: resources
      type: string
    - name: description
      type: string
    - name: files
      type:
        items: File
        type: array
    - name: reference__fasta__base
      type: File
    - name: config__algorithm__expression_caller
      type:
        items: string
        type: array
    - name: rgnames__lb
      type:
      - 'null'
      - string
    - name: rgnames__rg
      type: string
    - name: reference__hisat2__indexes
      type: File
    - name: config__algorithm__fusion_caller
      type:
        items: string
        type: array
    - name: config__algorithm__aligner
      type: string
    - name: rgnames__pl
      type: string
    - name: genome_build
      type: string
    - name: rgnames__pu
      type: string
    - name: genome_resources__rnaseq__transcripts
      type: File
    - name: config__algorithm__quality_format
      type: string
    - name: analysis
      type: string
    - name: rgnames__sample
      type: string
    - name: rgnames__lane
      type: string
    name: trim_rec
    type: record
outputs:
- id: fusion__fasta
  type: File
- id: fusion__json
  type: File
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
