$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=prep_rec:resources;description;files;reference__fasta__base;config__algorithm__expression_caller;rgnames__lb;rgnames__rg;reference__hisat2__indexes;config__algorithm__fusion_caller;config__algorithm__aligner;rgnames__pl;genome_build;rgnames__pu;genome_resources__rnaseq__transcripts;config__algorithm__quality_format;analysis;rgnames__sample;rgnames__lane
- sentinel_inputs=files:var,rgnames__sample:var,reference__fasta__base:var,genome_build:var,genome_resources__rnaseq__transcripts:var,analysis:var,rgnames__pl:var,rgnames__pu:var,rgnames__lane:var,rgnames__rg:var,rgnames__lb:var,reference__hisat2__indexes:var,config__algorithm__aligner:var,config__algorithm__expression_caller:var,config__algorithm__fusion_caller:var,config__algorithm__quality_format:var,resources:var,description:var
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- prepare_sample
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
  - package: picard
    specs:
    - https://anaconda.org/bioconda/picard
  - package: samtools
    specs:
    - https://anaconda.org/bioconda/samtools
  - package: pysam>
    specs:
    - https://anaconda.org/bioconda/pysam>
    version:
    - 0.13.0
inputs:
- id: files
  type:
    items: File
    type: array
- id: rgnames__sample
  type: string
- id: reference__fasta__base
  secondaryFiles:
  - .fai
  - ^.dict
  type: File
- id: genome_build
  type: string
- id: genome_resources__rnaseq__transcripts
  secondaryFiles:
  - .db
  type: File
- id: analysis
  type: string
- id: rgnames__pl
  type: string
- id: rgnames__pu
  type: string
- id: rgnames__lane
  type: string
- id: rgnames__rg
  type: string
- id: rgnames__lb
  type:
  - 'null'
  - string
- id: reference__hisat2__indexes
  type: File
- id: config__algorithm__aligner
  type: string
- id: config__algorithm__expression_caller
  type:
    items: string
    type: array
- id: config__algorithm__fusion_caller
  type:
    items: string
    type: array
- id: config__algorithm__quality_format
  type: string
- id: resources
  type: string
- id: description
  type: string
outputs:
- id: prep_rec
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
    name: prep_rec
    type: record
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
