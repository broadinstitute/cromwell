$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=multi-parallel
- sentinel_outputs=ensemble_rec:ensemble__vrn_file;ensemble__validate__summary;ensemble__batch_samples;ensemble__batch_id
- sentinel_inputs=ensemble_prep_rec:record
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- combine_calls
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 16
  outdirMin: 10240
  ramMin: 49152
  tmpdirMin: 0
- class: dx:InputResourceRequirement
  indirMin: 4687
- class: SoftwareRequirement
  packages:
  - package: bcbio-variation-recall
    specs:
    - https://anaconda.org/bioconda/bcbio-variation-recall
inputs:
- id: ensemble_prep_rec
  type:
    fields:
    - name: batch_id
      type: string
    - name: variants__calls
      type:
        items: File
        type: array
    - name: variants__variantcallers
      type:
        items: string
        type: array
    - name: resources
      type: string
    - name: description
      type: string
    - name: batch_samples
      type:
      - 'null'
      - items: string
        type: array
    - name: validate__summary
      type:
      - File
      - 'null'
    - name: validate__tp
      type:
      - File
      - 'null'
    - name: validate__fp
      type:
      - File
      - 'null'
    - name: validate__fn
      type:
      - File
      - 'null'
    - name: vrn_file
      type: File
    - name: reference__fasta__base
      type: File
    - name: config__algorithm__vcfanno
      type:
      - 'null'
      - items: 'null'
        type: array
    - name: config__algorithm__variantcaller
      type: string
    - name: config__algorithm__coverage_interval
      type:
      - string
      - 'null'
    - name: metadata__batch
      type: string
    - name: config__algorithm__min_allele_fraction
      type: double
    - name: reference__snpeff__GRCh37_75
      type: File
    - name: reference__genome_context
      type:
        items: File
        type: array
    - name: reference__rtg
      type: File
    - name: config__algorithm__validate
      type:
      - File
      - 'null'
    - name: config__algorithm__validate_regions
      type:
      - File
      - 'null'
    - name: genome_build
      type: string
    - name: metadata__phenotype
      type: string
    - name: genome_resources__aliases__human
      type:
      - string
      - 'null'
      - boolean
    - name: config__algorithm__tools_off
      type:
        items: string
        type: array
    - name: config__algorithm__ensemble
      type: string
    - name: analysis
      type: string
    - name: config__algorithm__tools_on
      type:
        items: string
        type: array
    - name: config__algorithm__effects
      type: string
    - name: config__algorithm__variant_regions
      type:
      - File
      - 'null'
    - name: genome_resources__aliases__ensembl
      type: string
    - name: config__algorithm__exclude_regions
      type:
      - 'null'
      - items: 'null'
        type: array
    - name: genome_resources__aliases__snpeff
      type: string
    - name: config__algorithm__variant_regions_merged
      type:
      - File
      - 'null'
    - name: regions__sample_callable
      type:
      - File
      - 'null'
    - name: config__algorithm__callable_regions
      type: File
    name: ensemble_prep_rec
    type: record
outputs:
- id: ensemble_rec
  type:
    fields:
    - name: ensemble__vrn_file
      type:
      - File
      - 'null'
    - name: ensemble__validate__summary
      type:
      - File
      - 'null'
    - name: ensemble__batch_samples
      type:
        items: string
        type: array
    - name: ensemble__batch_id
      type: string
    name: ensemble_rec
    type: record
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
