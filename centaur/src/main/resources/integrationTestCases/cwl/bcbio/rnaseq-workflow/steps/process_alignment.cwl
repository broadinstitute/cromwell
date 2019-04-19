$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=align_bam
- sentinel_inputs=trim_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- process_alignment
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 2
  outdirMin: 1028
  ramMin: 4096
  tmpdirMin: 2
- class: dx:InputResourceRequirement
  indirMin: 4
- class: SoftwareRequirement
  packages:
  - package: star
    specs:
    - https://anaconda.org/bioconda/star
  - package: hisat2
    specs:
    - https://anaconda.org/bioconda/hisat2
  - package: tophat
    specs:
    - https://anaconda.org/bioconda/tophat
  - package: samtools
    specs:
    - https://anaconda.org/bioconda/samtools
  - package: sambamba
    specs:
    - https://anaconda.org/bioconda/sambamba
  - package: seqtk
    specs:
    - https://anaconda.org/bioconda/seqtk
inputs:
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
- id: align_bam
  secondaryFiles:
  - .bai
  type: File
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
