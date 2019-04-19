$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=count_file,quant__tsv,quant__hdf5,quant__fusion
- sentinel_inputs=trim_rec:record,align_bam:var
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- rnaseq_quantitate
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-rnaseq:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 2
  outdirMin: 1025
  ramMin: 4096
  tmpdirMin: 1
- class: dx:InputResourceRequirement
  indirMin: 4
- class: SoftwareRequirement
  packages:
  - package: sailfish
    specs:
    - https://anaconda.org/bioconda/sailfish
  - package: salmon
    specs:
    - https://anaconda.org/bioconda/salmon
  - package: kallisto>
    specs:
    - https://anaconda.org/bioconda/kallisto>
    version:
    - 0.43.1
  - package: subread
    specs:
    - https://anaconda.org/bioconda/subread
  - package: gffread
    specs:
    - https://anaconda.org/bioconda/gffread
  - package: r
    specs:
    - https://anaconda.org/bioconda/r
    version:
    - 3.4.1
  - package: r-wasabi
    specs:
    - https://anaconda.org/bioconda/r-wasabi
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
- id: align_bam
  secondaryFiles:
  - .bai
  type: File
outputs:
- id: count_file
  type: File
- id: quant__tsv
  type: File
- id: quant__hdf5
  type: File
- id: quant__fusion
  type: File
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
