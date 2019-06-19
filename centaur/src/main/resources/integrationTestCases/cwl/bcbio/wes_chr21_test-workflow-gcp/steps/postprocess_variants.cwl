$namespaces:
  dx: https://www.dnanexus.com/cwl#
arguments:
- position: 0
  valueFrom: sentinel_runtime=cores,$(runtime['cores']),ram,$(runtime['ram'])
- sentinel_parallel=batch-single
- sentinel_outputs=vrn_file
- sentinel_inputs=batch_rec:record,vrn_file:var
- run_number=0
baseCommand:
- bcbio_nextgen.py
- runfn
- postprocess_variants
- cwl
class: CommandLineTool
cwlVersion: v1.0
hints:
- class: DockerRequirement
  dockerImageId: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
  dockerPull: quay.io/bcbio/bcbio-vc:1.1.4a-741877e
- class: ResourceRequirement
  coresMin: 16
  outdirMin: 10299
  ramMin: 49152
  tmpdirMin: 30
- class: dx:InputResourceRequirement
  indirMin: 24883
- class: SoftwareRequirement
  packages:
  - package: snpeff
    specs:
    - https://anaconda.org/bioconda/snpeff
    version:
    - 4.3.1t
inputs:
- id: batch_rec
  type:
    items:
      fields:
      - name: resources
        type: string
      - name: description
        type: string
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
      - name: genome_resources__variation__clinvar
        type: File
      - name: genome_resources__variation__esp
        type: File
      - name: metadata__batch
        type: string
      - name: genome_resources__variation__lcr
        type: File
      - name: config__algorithm__min_allele_fraction
        type: double
      - name: reference__snpeff__GRCh37_75
        type: File
      - name: vrn_file
        type:
        - 'null'
        - string
      - name: genome_resources__variation__train_hapmap
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
      - name: genome_resources__variation__1000g
        type: File
      - name: config__algorithm__validate_regions
        type:
        - File
        - 'null'
      - name: genome_build
        type: string
      - name: genome_resources__variation__exac
        type: File
      - name: genome_resources__variation__gnomad_exome
        type:
        - 'null'
        - string
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
      - name: genome_resources__variation__dbsnp
        type: File
      - name: genome_resources__variation__polyx
        type: File
      - name: genome_resources__variation__encode_blacklist
        type: File
      - name: genome_resources__variation__cosmic
        type: File
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
      - name: genome_resources__variation__train_indels
        type: File
      - name: genome_resources__aliases__snpeff
        type: string
      - name: align_bam
        type:
        - File
        - 'null'
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
      name: batch_rec
      type: record
    type: array
- id: vrn_file_toolinput
  secondaryFiles:
  - .tbi
  type: File
outputs:
- id: vrn_file
  secondaryFiles:
  - .tbi
  type: File
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(JSON.stringify(inputs))
    entryname: cwl.inputs.json
