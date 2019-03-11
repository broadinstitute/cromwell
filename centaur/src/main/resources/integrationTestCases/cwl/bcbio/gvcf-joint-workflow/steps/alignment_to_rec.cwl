$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=alignment_rec:resources;description;config__algorithm__align_split_size;files;config__algorithm__trim_reads;reference__fasta__base;config__algorithm__adapters;rgnames__lb;rgnames__rg;rgnames__lane;reference__bwa__indexes;config__algorithm__bam_clean;config__algorithm__aligner;rgnames__pl;rgnames__pu;config__algorithm__mark_duplicates;analysis;rgnames__sample;config__algorithm__variant_regions
- sentinel_inputs=files:var,analysis:var,config__algorithm__align_split_size:var,reference__fasta__base:var,rgnames__pl:var,rgnames__sample:var,rgnames__pu:var,rgnames__lane:var,rgnames__rg:var,rgnames__lb:var,reference__bwa__indexes:var,config__algorithm__aligner:var,config__algorithm__trim_reads:var,config__algorithm__adapters:var,config__algorithm__bam_clean:var,config__algorithm__variant_regions:var,config__algorithm__mark_duplicates:var,resources:var,description:var
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- alignment_to_rec
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 1
  outdirMin: 1028
  ramMin: 2048
  tmpdirMin: 2
- class: dx:InputResourceRequirement
  indirMin: 0
inputs:
- id: files
  type:
    items:
      items: File
      type: array
    type: array
- id: analysis
  type:
    items: string
    type: array
- id: config__algorithm__align_split_size
  type:
    items:
    - 'null'
    - string
    type: array
- id: reference__fasta__base
  secondaryFiles:
  - .fai
  - ^.dict
  type:
    items: File
    type: array
- id: rgnames__pl
  type:
    items: string
    type: array
- id: rgnames__sample
  type:
    items: string
    type: array
- id: rgnames__pu
  type:
    items: string
    type: array
- id: rgnames__lane
  type:
    items: string
    type: array
- id: rgnames__rg
  type:
    items: string
    type: array
- id: rgnames__lb
  type:
    items:
    - 'null'
    - string
    type: array
- id: reference__bwa__indexes
  secondaryFiles:
  - ^.ann
  - ^.pac
  - ^.sa
  - ^.bwt
  type:
    items: File
    type: array
- id: config__algorithm__aligner
  type:
    items: string
    type: array
- id: config__algorithm__trim_reads
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: config__algorithm__adapters
  type:
    items:
      items: string
      type: array
    type: array
- id: config__algorithm__bam_clean
  type:
    items:
    - string
    - 'null'
    - boolean
    type: array
- id: config__algorithm__variant_regions
  type:
    items: File
    type: array
- id: config__algorithm__mark_duplicates
  type:
    items:
    - string
    - 'null'
    - boolean
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
- id: alignment_rec
  type:
    items:
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
          items: string
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
    type: array
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
