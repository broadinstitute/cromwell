$namespaces:
  arv: http://arvados.org/cwl#
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-combined
- sentinel_outputs=config__algorithm__callable_regions,config__algorithm__non_callable_regions,config__algorithm__callable_count
- sentinel_inputs=regions__callable:var,regions__nblock:var,metadata__batch:var,config__algorithm__nomap_split_size:var,config__algorithm__nomap_split_targets:var,reference__fasta__base:var,resources:var,description:var
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- combine_sample_regions
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 1
  outdirMin: 1025
  ramMin: 3072
  tmpdirMin: 1
- class: dx:InputResourceRequirement
  indirMin: 1
- class: SoftwareRequirement
  packages:
  - package: bedtools
    specs:
    - https://anaconda.org/bioconda/bedtools
  - package: htslib
    specs:
    - https://anaconda.org/bioconda/htslib
  - package: gatk4
    specs:
    - https://anaconda.org/bioconda/gatk4
- class: arv:APIRequirement
inputs:
- id: regions__callable
  type:
    items:
    - File
    - 'null'
    type: array
- id: regions__nblock
  type:
    items:
    - File
    - 'null'
    type: array
- id: metadata__batch
  type:
    items: string
    type: array
- id: config__algorithm__nomap_split_size
  type:
    items: long
    type: array
- id: config__algorithm__nomap_split_targets
  type:
    items: long
    type: array
- id: reference__fasta__base
  secondaryFiles:
  - .fai
  - ^.dict
  type:
    items: File
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
- id: config__algorithm__callable_regions
  type:
    items: File
    type: array
- id: config__algorithm__non_callable_regions
  type:
    items: File
    type: array
- id: config__algorithm__callable_count
  type:
    items: int
    type: array
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
